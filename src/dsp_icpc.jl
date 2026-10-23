# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    dsp_icpc(data::Q, config::DSPConfig, τ::Quantity{T}, pars_filter::PropDict) where {Q <: Table, T<:Real}

DSP for ICPC detectors. It needs the decay time `τ` of the detector and the filter parameters `pars_filter` for the optimal filter parameters
for the Trap, CUSP and ZAC filter. 
# Input data
The input data is a table with the following columns:
- `waveform`: waveform data
- `baseline`: baseline from FADC
- `timestamp`: timestamp from FADC
- `eventnumber`: event ID from FADC
- `daqenergy`: energy from FADC

# Output data
The output data is a table with the following columns:
- `blmean`: baseline mean
- `blsigma`: baseline sigma
- `blslope`: baseline slope
- `bloffset`: baseline offset
- `tailmean`: tail mean after PZ correction
- `tailsigma`: tail sigma after PZ correction
- `tailslope`: tail slope after PZ correction
- `tailoffset`: tail offset after PZ correction
- `t0`: start time of waveform drift
- `t10`: timepoint of 10% of the pulse amplitude (maximum of the 101010 trapezoid)
- `t50`: timepoint of 50% of the pulse amplitude (maximum of the 101010 trapezoid)
- `t80`: timepoint of 80% of the pulse amplitude (maximum of the 101010 trapezoid)
- `t90`: timepoint of 90% of the pulse amplitude (maximum of the 101010 trapezoid)
- `t99`: timepoint of 99% of the pulse amplitude (maximum of the 101010 trapezoid)
- `t50_current`: timepoint of current rise to 50% of maximum
- `t_amax`: quadratically interpolated timepoint of maximum current
- `tail_τ`: tail decay time
- `tail_mean`: tail mean before PZ correction
- `tail_sigma`: tail sigma before PZ correction
- `e_max`: maximum of waveform
- `e_min`: minimum of waveform
- `e_10410`: energy of waveform with trapezoidal filter of 10µs rise time with 4µs flat-top
- `e_313`: energy of waveform with trapezoidal filter of 3µs rise time with 1µs flat-top
- `e_101010`: energy of waveform with trapezoidal filter of 10µs rise time with 10µs flat-top
- `a_sg`: maximum of the current from the Savitzky-Golay derivative (`sg_flt_degree`) with the optimized window length from the `aoeopt` pars
- `a_sg_<ns>`: the same with the fixed window lengths `a_sg_windows` of the config
- `a_raw`: maximum of the raw current (simple derivative)
- `a_mwa`: maximum of the raw current after the multi-pass moving-window average (`a_mwa_passes`, alternating left and right) with the optimized window length from the `aoeopt` pars
- `a_mwa_<ns>`: the same with the fixed window lengths `a_mwa_windows` of the config
- `t_a_sg`, `t_a_mwa`: quadratically interpolated timepoints of the `a_sg` and `a_mwa` maxima
- `e_10410_inv`: maximum of inverted waveform with trapezoidal filter of 10µs rise time with 4µs flat-top
- `e_313_inv`: maximum of inverted waveform with trapezoidal filter of 3µs rise time with 1µs flat-top
- `t0_inv`: start time of inverted waveform drift
- `e_trap`: energy of waveform with trapezoidal filter of optimized rise and flat-top time
- `e_cusp`: energy of waveform with CUSP filter of optimized rise and flat-top time
- `e_zac`: energy of waveform with ZAC filter of optimized rise and flat-top time
- `qdrift`: Q-drift parameter
- `lq80_<early><gap><late>`, `lq90_…`, `lqamax_…`: LQ parameters from `t80`, `t90` and `t_amax`, charge in the late window minus charge in the early window, one column per `[early, gap, late]` of `lq_windows` in the config (lengths in 0.1 µs, e.g. `lq80_250025`: 2.5 µs early, no gap, 2.5 µs late)
- `blfc`: baseline from FADC
- `timestamp`: timestamp from FADC
- `eventID_fadc`: event ID from FADC
- `e_fc`: energy from FADC
- `pretrace_diff`: difference between first sample and baseline mean
- `drift_time`: drift time between t0 and 90% of waveform maximum
- `inTrace_intersect`: position of in-trace pile-up
- `inTrace_n`: multiplicity of in-trace pile-up
- `n_sat_low`: number of samples the waveform is saturated at low of FADC range
- `n_sat_high`: number of samples the waveform is saturated at high of FADC range
- `n_sat_low_cons`: number of consecutive samples the waveform is saturated at low of FADC range
- `n_sat_high_cons`: number of consecutive samples the waveform is saturated at high of FADC range
"""
function dsp_icpc(data::Q, config::DSPConfig, τ::Quantity{T}, pars_filter::PropDict; f_evaluate_qc::Union{Function, Missing}=missing) where {Q <: Table, T<:Real}
    # get config parameters
    bl_window                = config.bl_window
    t0_threshold             = config.t0_threshold
    tail_window              = config.tail_window
    inTraceCut_std_threshold = config.inTraceCut_std_threshold
    sg_flt_degree            = config.sg_flt_degree
    current_window           = config.current_window
    qdrift_int_length        = config.qdrift_int_length
    
    # get optimal filter parameters
    trap_rt, trap_ft = get_fltpars(pars_filter, :trap, config)
    cusp_rt, cusp_ft = get_fltpars(pars_filter, :cusp, config)
    zac_rt, zac_ft = get_fltpars(pars_filter, :zac, config)
    sg_wl   = get_fltpars(pars_filter, :sg, config)
    mwa_wl  = get_fltpars(pars_filter, :mwa, config)

    # get waveform data
    data = decode_data(data)
    wvfs = data.waveform
    blfc = data.baseline
    ts   = data.timestamp
    evID = data.eventnumber
    efc  = data.daqenergy

    # get CUSP and ZAC filter length and flt scale
    flt_length_zac              = config.flt_length_zac
    zac_scale                   = ustrip(NoUnits, flt_length_zac/step(wvfs[1].time))
    flt_length_cusp             = config.flt_length_cusp
    cusp_scale                  = ustrip(NoUnits, flt_length_cusp/step(wvfs[1].time))

    # get number of samples the waveform is saturated at low and high of FADC range
    bit_depth = config.kwargs_pars.fc_bit_depth # of FlashCam FADC
    sat_low, sat_high = 0, 2^bit_depth - bit_depth
    sat_stats = saturation.(wvfs, sat_low, sat_high)

    # set τ for CUSP filter to very high number to switch of CR filter
    τ_cusp = 10000000.0u"µs"
    τ_zac = 10000000.0u"µs"

    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs, leftendpoint(bl_window), rightendpoint(bl_window))

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)

    # get QC classifier labels (skip expensive ML classification if no model)
    qc_labels = !ismissing(f_evaluate_qc) ? Int.(get_qc_classifier(wvfs, f_evaluate_qc)) : fill(-1, length(wvfs))
    
    # get wvf maximum
    wvf_max = maximum.(wvfs.signal)
    wvf_min = minimum.(wvfs.signal)

    # extract decay times
    tail_stats = tailstats.(wvfs, leftendpoint(tail_window), rightendpoint(tail_window))

    # deconvolute waveform 
    # --> wvfs = wvfs_pz
    deconv_flt = InvCRFilter(τ)
    wvfs = deconv_flt.(wvfs)

    # get tail mean, std and slope
    pz_stats = signalstats.(wvfs, leftendpoint(tail_window), rightendpoint(tail_window))

    # t0 determination
    t0 = get_t0(wvfs, t0_threshold; flt_pars=config.kwargs_pars.t0_flt_pars, mintot=config.kwargs_pars.t0_mintot)

    # if all waveforms are saturated set threshold to 1.0 to avoid numerical problems
    # replace!(wvf_max, zero(wvf_max[1]) => one(wvf_max[1]))

    # robust energy reconstruction with long, middle and short rise and flat-top times
    uflt_10410 = TrapezoidalChargeFilter(10u"µs", 4u"µs")
    e_10410  = maximum.((uflt_10410.(wvfs)).signal)

    uflt_535 = TrapezoidalChargeFilter(5u"µs", 3u"µs")
    e_535  = maximum.((uflt_535.(wvfs)).signal)

    uflt_313 = TrapezoidalChargeFilter(3u"µs", 1u"µs")
    e_313  = maximum.((uflt_313.(wvfs)).signal)

    uflt_101010 = TrapezoidalChargeFilter(10u"µs", 10u"µs")
    e_101010  = maximum.((uflt_101010.(wvfs)).signal)

    # get threshold points in rise
    t10 = get_threshold(wvfs, e_101010 .* 0.1; mintot=config.kwargs_pars.tx_mintot)
    t50 = get_threshold(wvfs, e_101010 .* 0.5; mintot=config.kwargs_pars.tx_mintot)
    t80 = get_threshold(wvfs, e_101010 .* 0.8; mintot=config.kwargs_pars.tx_mintot)
    t90 = get_threshold(wvfs, e_101010 .* 0.9; mintot=config.kwargs_pars.tx_mintot)
    t99 = get_threshold(wvfs, e_101010 .* 0.99; mintot=config.kwargs_pars.tx_mintot)
    
    drift_time = uconvert.(u"ns", t90 - t0)

    # get t_amax with quadratic interpolation around the maximum current sample
    wvfs_deriv = DerivativeFilter(1).(wvfs)
    t_amax = get_wvf_maximum.(wvfs_deriv, first(wvfs_deriv[1].time), last(wvfs_deriv[1].time)).t

    # get Q-drift parameter
    qdrift = get_qdrift(wvfs, t0, qdrift_int_length; pol_power=config.kwargs_pars.int_interpolation_order, sign_est_length=config.kwargs_pars.int_interpolation_length)

    # get LQ parameters
    wvfs_int = IntegratorFilter(1).(wvfs)
    lq80   = NamedTuple(_lq_column_name(:t80, w)    => get_lq(wvfs_int, t80, w...; pol_power=config.kwargs_pars.int_interpolation_order, sign_est_length=config.kwargs_pars.int_interpolation_length) for w in config.lq_windows)
    lq90   = NamedTuple(_lq_column_name(:t90, w)    => get_lq(wvfs_int, t90, w...; pol_power=config.kwargs_pars.int_interpolation_order, sign_est_length=config.kwargs_pars.int_interpolation_length) for w in config.lq_windows)
    lqamax = NamedTuple(_lq_column_name(:t_amax, w) => get_lq(wvfs_int, t_amax, w...; pol_power=config.kwargs_pars.int_interpolation_order, sign_est_length=config.kwargs_pars.int_interpolation_length) for w in config.lq_windows)

    # signal estimator for precise energy reconstruction
    signal_estimator = SignalEstimator(PolynomialDNI(config.kwargs_pars.sig_interpolation_order, config.kwargs_pars.sig_interpolation_length))

    # get trap energy of optimized rise and flat-top time
    uflt_trap_rtft = TrapezoidalChargeFilter(trap_rt, trap_ft)
    wvfs_flt = uflt_trap_rtft.(wvfs)

    e_trap = signal_estimator.(wvfs_flt, t50 .+ (trap_rt + trap_ft/2))
    e_trap_extremestats = extremestats.(wvfs_flt)

    # get cusp energy of optimized rise and flat-top time
    uflt_cusp_rtft = CUSPChargeFilter(cusp_rt, cusp_ft, τ_cusp, flt_length_cusp, cusp_scale)
    wvfs_flt = uflt_cusp_rtft.(wvfs)

    e_cusp = signal_estimator.(wvfs_flt, t50 .+ (flt_length_cusp /2))
    e_cusp_extremestats = extremestats.(wvfs_flt)

    # get zac energy of optimized rise and flat-top time
    uflt_zac_rtft = ZACChargeFilter(zac_rt, zac_ft, τ_zac, flt_length_zac, zac_scale)
    wvfs_flt = uflt_zac_rtft.(wvfs)
    
    e_zac = signal_estimator.(uflt_zac_rtft.(wvfs), t50 .+ (flt_length_zac /2))
    e_zac_extremestats = extremestats.(wvfs_flt)

    # extract current with optimal SG filter length with second order polynominal and first derivative
    wvfs_sgflt_deriv = SavitzkyGolayFilter(sg_wl, sg_flt_degree, 1).(wvfs)
    a_sg_max = get_a(wvfs_sgflt_deriv, current_window)
    a_sg, t_a_sg = a_sg_max.max, a_sg_max.t

    a_raw = get_a(wvfs_deriv, current_window).max
    a_mwa_max = get_a(MovingWindowMultiFilter(mwa_wl; passes = config.a_mwa_passes).(wvfs_deriv), current_window)
    a_mwa, t_a_mwa = a_mwa_max.max, a_mwa_max.t

    a_sg_fixed  = NamedTuple(Symbol(:a_sg_,  round(Int, ustrip(u"ns", wl))) => get_a(SavitzkyGolayFilter(wl, sg_flt_degree, 1).(wvfs), current_window).max for wl in config.a_sg_windows)
    a_mwa_fixed = NamedTuple(Symbol(:a_mwa_, round(Int, ustrip(u"ns", wl))) => get_a(MovingWindowMultiFilter(wl; passes = config.a_mwa_passes).(wvfs_deriv), current_window).max for wl in config.a_mwa_windows)

    # get in-trace pile-up
    inTrace_pileUp = get_triggers(wvfs_sgflt_deriv, inTraceCut_std_threshold, bl_window; mintot=config.kwargs_pars.intrace_mintot)
    
    # get position of current rise
    thres = maximum.(wvfs_sgflt_deriv.signal) .* 0.5
    # replace!(thres, zero(thres[1]) => one(thres[1]))

    t50_current = get_threshold(wvfs_sgflt_deriv, thres; mintot=config.kwargs_pars.tx_mintot)

    # invert waveform for DC tagging
    # wvfs --> wvfs_pz_inv
    wvfs = multiply_waveform.(wvfs, -1.0)

    # get inverted waveform maximum for long and short filter
    e_10410_max_inv  = maximum.(uflt_10410.(wvfs).signal)

    e_313_max_inv  = maximum.(uflt_313.(wvfs).signal)

    # t0 determination
    t0_inv = get_t0(wvfs, t0_threshold; mintot=config.kwargs_pars.t0_mintot)

    # output Table 
    return TypedTables.Table(; blmean = bl_stats.mean, blsigma = bl_stats.sigma, blslope = bl_stats.slope, bloffset = bl_stats.offset, 
    tailmean = pz_stats.mean, tailsigma = pz_stats.sigma, tailslope = pz_stats.slope, tailoffset = pz_stats.offset,
    qc_label = qc_labels,
    t0 = t0, t10 = t10, t50 = t50, t80 = t80, t90 = t90, t99 = t99,
    t50_current = t50_current, t_amax = t_amax, t_a_sg = t_a_sg, t_a_mwa = t_a_mwa,
    drift_time = drift_time,
    tail_τ = tail_stats.τ, tail_mean = tail_stats.mean, tail_sigma = tail_stats.sigma,
    e_max = wvf_max, e_min = wvf_min,
    e_10410 = e_10410, e_535 = e_535, e_313 = e_313, e_101010 = e_101010,
    e_10410_inv = e_10410_max_inv, e_313_inv = e_313_max_inv,
    t0_inv = t0_inv,
    e_trap = e_trap, e_cusp = e_cusp, e_zac = e_zac,
    e_trap_max = e_trap_extremestats.max, e_cusp_max = e_cusp_extremestats.max, e_zac_max = e_zac_extremestats.max,
    t_trap_max = e_trap_extremestats.tmax, t_cusp_max = e_cusp_extremestats.tmax, t_zac_max = e_zac_extremestats.tmax,
    qdrift = qdrift, lq80..., lq90..., lqamax...,
    a_sg = a_sg, a_raw = a_raw, a_mwa = a_mwa, a_sg_fixed..., a_mwa_fixed...,
    blfc = blfc, timestamp = ts, eventID_fadc = evID, e_fc = efc,
    inTrace_intersect = inTrace_pileUp.intersect, inTrace_n = inTrace_pileUp.n,
    n_sat_low = sat_stats.low, n_sat_high = sat_stats.high, n_sat_low_cons = sat_stats.max_cons_low, n_sat_high_cons = sat_stats.max_cons_high
    )
end
export dsp_icpc


"""
    dsp_icpc_compressed(data::Q, config::DSPConfig, τ::Quantity{T}, pars_filter::PropDict) where {Q <: Table, T<:Real}

DSP for ICPC detectors. It needs the decay time `τ` of the detector and the filter parameters `pars_filter` for the optimal filter parameters
for the Trap, CUSP and ZAC filter. 
# Input data
The input data is a table with the following columns:
- `waveform`: waveform data
- `baseline`: baseline from FADC
- `timestamp`: timestamp from FADC
- `eventnumber`: event ID from FADC
- `daqenergy`: energy from FADC

# Output data
The output data is a table with the following columns:
- `blmean`: baseline mean
- `blsigma`: baseline sigma
- `blslope`: baseline slope
- `bloffset`: baseline offset
- `tailmean`: tail mean after PZ correction
- `tailsigma`: tail sigma after PZ correction
- `tailslope`: tail slope after PZ correction
- `tailoffset`: tail offset after PZ correction
- `t0`: start time of waveform drift
- `t10`: timepoint of 10% of the pulse amplitude (maximum of the 101010 trapezoid)
- `t50`: timepoint of 50% of the pulse amplitude (maximum of the 101010 trapezoid)
- `t80`: timepoint of 80% of the pulse amplitude (maximum of the 101010 trapezoid)
- `t90`: timepoint of 90% of the pulse amplitude (maximum of the 101010 trapezoid)
- `t99`: timepoint of 99% of the pulse amplitude (maximum of the 101010 trapezoid)
- `t50_current`: timepoint of current rise to 50% of maximum
- `t_amax`: quadratically interpolated timepoint of maximum current
- `tail_τ`: tail decay time
- `tail_mean`: tail mean before PZ correction
- `tail_sigma`: tail sigma before PZ correction
- `e_max`: maximum of waveform
- `e_min`: minimum of waveform
- `e_10410`: energy of waveform with trapezoidal filter of 10µs rise time with 4µs flat-top
- `e_313`: energy of waveform with trapezoidal filter of 3µs rise time with 1µs flat-top
- `e_101010`: energy of waveform with trapezoidal filter of 10µs rise time with 10µs flat-top
- `a_sg`: maximum of the current from the Savitzky-Golay derivative (`sg_flt_degree`) with the optimized window length from the `aoeopt` pars
- `a_sg_<ns>`: the same with the fixed window lengths `a_sg_windows` of the config
- `a_raw`: maximum of the raw current (simple derivative)
- `a_mwa`: maximum of the raw current after the multi-pass moving-window average (`a_mwa_passes`, alternating left and right) with the optimized window length from the `aoeopt` pars
- `a_mwa_<ns>`: the same with the fixed window lengths `a_mwa_windows` of the config
- `t_a_sg`, `t_a_mwa`: quadratically interpolated timepoints of the `a_sg` and `a_mwa` maxima
- `e_10410_inv`: maximum of inverted waveform with trapezoidal filter of 10µs rise time with 4µs flat-top
- `e_313_inv`: maximum of inverted waveform with trapezoidal filter of 3µs rise time with 1µs flat-top
- `t0_inv`: start time of inverted waveform drift
- `e_trap`: energy of waveform with trapezoidal filter of optimized rise and flat-top time
- `e_cusp`: energy of waveform with CUSP filter of optimized rise and flat-top time
- `e_zac`: energy of waveform with ZAC filter of optimized rise and flat-top time
- `qdrift`: Q-drift parameter
- `lq80_<early><gap><late>`, `lq90_…`, `lqamax_…`: LQ parameters from `t80`, `t90` and `t_amax`, charge in the late window minus charge in the early window, one column per `[early, gap, late]` of `lq_windows` in the config (lengths in 0.1 µs, e.g. `lq80_250025`: 2.5 µs early, no gap, 2.5 µs late)
- `blfc`: baseline from FADC
- `timestamp`: timestamp from FADC
- `eventID_fadc`: event ID from FADC
- `e_fc`: energy from FADC
- `pretrace_diff`: difference between first sample and baseline mean
- `drift_time`: drift time between t0 and 90% of waveform maximum
- `inTrace_intersect`: position of in-trace pile-up
- `inTrace_n`: multiplicity of in-trace pile-up
- `n_sat_low`: number of samples the waveform is saturated at low of FADC range
- `n_sat_high`: number of samples the waveform is saturated at high of FADC range
- `n_sat_low_cons`: number of consecutive samples the waveform is saturated at low of FADC range
- `n_sat_high_cons`: number of consecutive samples the waveform is saturated at high of FADC range
"""
function dsp_icpc_compressed(data::Q, config::DSPConfig, τ::Quantity{T}, pars_filter::PropDict; f_evaluate_qc::Union{Function, Missing}=missing) where {Q <: Table, T<:Real}
    # get config parameters
    bl_window                = config.bl_window
    t0_threshold             = config.t0_threshold
    tail_window              = config.tail_window
    inTraceCut_std_threshold = config.inTraceCut_std_threshold
    sg_flt_degree            = config.sg_flt_degree
    current_window           = config.current_window
    qdrift_int_length        = config.qdrift_int_length
    
    # get optimal filter parameters
    trap_rt, trap_ft = get_fltpars(pars_filter, :trap, config)
    cusp_rt, cusp_ft = get_fltpars(pars_filter, :cusp, config)
    zac_rt, zac_ft = get_fltpars(pars_filter, :zac, config)
    sg_wl   = get_fltpars(pars_filter, :sg, config)
    mwa_wl  = get_fltpars(pars_filter, :mwa, config)

    # get waveform data
    data = decode_data(data)
    wvfs_pre = data.waveform_presummed
    wvfs_wdw = data.waveform_windowed
    presum_rate = data.presum_rate
    blfc = data.baseline
    ts   = data.timestamp
    evID = data.eventnumber
    efc  = data.daqenergy
    t_sat_lo = data.t_sat_lo
    t_sat_hi = data.t_sat_hi
    deadtime = data.deadtime

    presum_rate_value = only(unique(presum_rate))

    # get CUSP and ZAC filter length and flt scale
    flt_length_zac              = config.flt_length_zac
    zac_scale                   = ustrip(NoUnits, flt_length_zac/step(wvfs_pre[1].time))
    flt_length_cusp             = config.flt_length_cusp
    cusp_scale                  = ustrip(NoUnits, flt_length_cusp/step(wvfs_pre[1].time))

    # get number of samples the waveform is saturated at low and high of FADC range
    bit_depth = config.kwargs_pars.fc_bit_depth # of FlashCam FADC
    sat_low, sat_high = 0, (2^bit_depth - bit_depth) * first(presum_rate_value)
    sat_stats = saturation.(wvfs_pre, sat_low, sat_high)

    # set τ for CUSP filter to very high number to switch of CR filter
    τ_cusp = 10000000.0u"µs"
    τ_zac = 10000000.0u"µs"

    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs_pre, leftendpoint(bl_window), rightendpoint(bl_window))

    # substract baseline from waveforms
    wvfs_pre = shift_waveform.(wvfs_pre, -bl_stats.mean)
    wvfs_wdw = shift_waveform.(wvfs_wdw, -bl_stats.mean ./ presum_rate_value)

    # get QC classifier labels (skip expensive ML classification if no model)
    qc_labels = !ismissing(f_evaluate_qc) ? Int.(get_qc_classifier_compressed(wvfs_pre, f_evaluate_qc)) : fill(-1, length(wvfs_pre))
    
    # get wvf extremestats
    wvf_pre_extremestats = extremestats.(wvfs_pre)
    wvf_wdw_extremestats = extremestats.(wvfs_wdw)

    # extract decay times
    tail_stats = tailstats.(wvfs_pre, leftendpoint(tail_window), rightendpoint(tail_window))

    # deconvolute waveform 
    # --> wvfs = wvfs_pz
    deconv_flt = InvCRFilter(τ)
    wvfs_pre = deconv_flt.(wvfs_pre)
    wvfs_wdw = deconv_flt.(wvfs_wdw)

    # auxiliary pz-corrected waveform window determination
    auxpz1_stats = signalstats.(wvfs_pre, leftendpoint(config.auxpz1_window), rightendpoint(config.auxpz1_window))
    auxpz2_stats = signalstats.(wvfs_pre, leftendpoint(config.auxpz2_window), rightendpoint(config.auxpz2_window))

    # get tail mean, std and slope
    pz_stats = signalstats.(wvfs_pre, leftendpoint(tail_window), rightendpoint(tail_window))

    # t0 determination
    t0 = get_t0(wvfs_wdw, t0_threshold; flt_pars=config.kwargs_pars.t0_flt_pars, mintot=config.kwargs_pars.t0_mintot)

    # robust energy reconstruction with long, middle and short rise and flat-top times
    uflt_10410 = TrapezoidalChargeFilter(10u"µs", 4u"µs")
    e_10410_stats  = extremestats.(uflt_10410.(wvfs_pre))
    
    uflt_535 = TrapezoidalChargeFilter(5u"µs", 3u"µs")
    e_535_stats  = extremestats.(uflt_535.(wvfs_pre))
    
    uflt_313 = TrapezoidalChargeFilter(3u"µs", 1u"µs")
    e_313_stats  = extremestats.(uflt_313.(wvfs_pre))

    uflt_101010 = TrapezoidalChargeFilter(10u"µs", 10u"µs")
    e_101010_stats  = extremestats.(uflt_101010.(wvfs_pre))

    # get threshold points in rise
    wvf_amp_wdw = e_101010_stats.max ./ presum_rate_value
    t10 = get_threshold(wvfs_wdw, wvf_amp_wdw .* 0.1; mintot=config.kwargs_pars.tx_mintot)
    t50 = get_threshold(wvfs_wdw, wvf_amp_wdw .* 0.5; mintot=config.kwargs_pars.tx_mintot)
    t50_pre = get_threshold(wvfs_pre, e_101010_stats.max .* 0.5; mintot=config.kwargs_pars.tx_mintot)
    t80 = get_threshold(wvfs_wdw, wvf_amp_wdw .* 0.8; mintot=config.kwargs_pars.tx_mintot)
    t90 = get_threshold(wvfs_wdw, wvf_amp_wdw .* 0.9; mintot=config.kwargs_pars.tx_mintot)
    t99 = get_threshold(wvfs_wdw, wvf_amp_wdw .* 0.99; mintot=config.kwargs_pars.tx_mintot)
    
    drift_time = uconvert.(u"ns", t90 - t0)

    # get t_amax with quadratic interpolation around the maximum current sample
    wvfs_deriv = DerivativeFilter(1).(wvfs_wdw)
    t_amax = get_wvf_maximum.(wvfs_deriv, first(wvfs_deriv[1].time), last(wvfs_deriv[1].time)).t

    # get Q-drift parameter
    qdrift = get_qdrift(wvfs_wdw, t0, qdrift_int_length; pol_power=config.kwargs_pars.int_interpolation_order, sign_est_length=config.kwargs_pars.int_interpolation_length)

    # get LQ parameters
    wvfs_int = IntegratorFilter(1).(wvfs_wdw)
    lq80   = NamedTuple(_lq_column_name(:t80, w)    => get_lq(wvfs_int, t80, w...; pol_power=config.kwargs_pars.int_interpolation_order, sign_est_length=config.kwargs_pars.int_interpolation_length) for w in config.lq_windows)
    lq90   = NamedTuple(_lq_column_name(:t90, w)    => get_lq(wvfs_int, t90, w...; pol_power=config.kwargs_pars.int_interpolation_order, sign_est_length=config.kwargs_pars.int_interpolation_length) for w in config.lq_windows)
    lqamax = NamedTuple(_lq_column_name(:t_amax, w) => get_lq(wvfs_int, t_amax, w...; pol_power=config.kwargs_pars.int_interpolation_order, sign_est_length=config.kwargs_pars.int_interpolation_length) for w in config.lq_windows)

    # signal estimator for precise energy reconstruction
    signal_estimator = SignalEstimator(PolynomialDNI(config.kwargs_pars.sig_interpolation_order, config.kwargs_pars.sig_interpolation_length))

    # get trap energy of optimized rise and flat-top time
    uflt_trap_rtft = TrapezoidalChargeFilter(trap_rt, trap_ft)
    wvfs_flt = uflt_trap_rtft.(wvfs_pre)

    e_trap = signal_estimator.(wvfs_flt, t50_pre .+ (trap_rt + trap_ft/2))
    e_trap_extremestats = extremestats.(wvfs_flt)

    # get cusp energy of optimized rise and flat-top time
    uflt_cusp_rtft = CUSPChargeFilter(cusp_rt, cusp_ft, τ_cusp, flt_length_cusp, cusp_scale)
    wvfs_flt = uflt_cusp_rtft.(wvfs_pre)

    e_cusp = signal_estimator.(wvfs_flt, t50_pre .+ (flt_length_cusp /2))
    e_cusp_extremestats = extremestats.(wvfs_flt)

    # get zac energy of optimized rise and flat-top time
    uflt_zac_rtft = ZACChargeFilter(zac_rt, zac_ft, τ_zac, flt_length_zac, zac_scale)
    wvfs_flt = uflt_zac_rtft.(wvfs_pre)

    e_zac = signal_estimator.(wvfs_flt, t50_pre .+ (flt_length_zac /2))
    e_zac_extremestats = extremestats.(wvfs_flt)

    # extract current with optimal SG filter length with second order polynominal and first derivative
    a_raw = get_a(wvfs_deriv, current_window).max
    a_mwa_max = get_a(MovingWindowMultiFilter(mwa_wl; passes = config.a_mwa_passes).(wvfs_deriv), current_window)
    a_sg_max  = get_a(SavitzkyGolayFilter(sg_wl, sg_flt_degree, 1).(wvfs_wdw), current_window)
    a_mwa, t_a_mwa = a_mwa_max.max, a_mwa_max.t
    a_sg,  t_a_sg  = a_sg_max.max,  a_sg_max.t

    a_sg_fixed  = NamedTuple(Symbol(:a_sg_,  round(Int, ustrip(u"ns", wl))) => get_a(SavitzkyGolayFilter(wl, sg_flt_degree, 1).(wvfs_wdw), current_window).max for wl in config.a_sg_windows)
    a_mwa_fixed = NamedTuple(Symbol(:a_mwa_, round(Int, ustrip(u"ns", wl))) => get_a(MovingWindowMultiFilter(wl; passes = config.a_mwa_passes).(wvfs_deriv), current_window).max for wl in config.a_mwa_windows)

    # Trigger filter with trapezoidal filter
    trig_flt = TrapezoidalChargeFilter(2u"µs", 1u"µs")
    wvfs_trigflt = trig_flt.(wvfs_pre)
    triggers = get_triggers(wvfs_trigflt, inTraceCut_std_threshold, bl_window; mintot=2.0u"µs")

    # get position of current rise
    wvfs_sgflt_deriv = SavitzkyGolayFilter(sg_wl * presum_rate_value / 2, sg_flt_degree, 1).(wvfs_pre)
    thres = maximum.(wvfs_sgflt_deriv.signal) .* 0.5

    t50_current = get_threshold(wvfs_sgflt_deriv, thres; mintot=config.kwargs_pars.tx_mintot)

    # invert waveform for DC tagging
    # wvfs --> wvfs_pz_inv
    wvfs_pre = multiply_waveform.(wvfs_pre, -1.0)
    wvfs_wdw = multiply_waveform.(wvfs_wdw, -1.0)

    # get inverted waveform maximum for long and short filter
    e_10410_max_inv  = maximum.(uflt_10410.(wvfs_pre).signal)

    e_313_max_inv  = maximum.(uflt_313.(wvfs_pre).signal)

    # inv t0 determination
    t0_inv = get_t0(wvfs_wdw, t0_threshold; mintot=config.kwargs_pars.t0_mintot)

    # output Table 
    return TypedTables.Table(;
        # fadc parameters
        blfc=blfc, timestamp=ts, eventID_fadc=evID, e_fc=efc, deadtime=deadtime,
        # saturation parameters 
        n_sat_low=sat_stats.low, n_sat_high=sat_stats.high, n_sat_low_cons=sat_stats.max_cons_low, n_sat_high_cons=sat_stats.max_cons_high,
        t_sat_lo=t_sat_lo, t_sat_hi=t_sat_hi,
        # baseline 
        blmean = bl_stats.mean, blsigma = bl_stats.sigma, blslope = bl_stats.slope, bloffset = bl_stats.offset, bl_slope_sigma = bl_stats.slope_residual_sigma,
        # ML QC label
        qc_label = qc_labels,
        # waveform extrema and timepoints
        e_max = wvf_pre_extremestats.max, e_min = wvf_wdw_extremestats.min, t_max = wvf_pre_extremestats.tmax, t_min = wvf_wdw_extremestats.tmin,
        e_max_pre = wvf_pre_extremestats.max, e_min_pre = wvf_pre_extremestats.min, t_max_pre = wvf_pre_extremestats.tmax, t_min_pre = wvf_pre_extremestats.tmin,
        # tail parameters
        tailmean = pz_stats.mean, tailsigma = pz_stats.sigma, tailslope = pz_stats.slope, tailoffset = pz_stats.offset,
        tail_τ = tail_stats.τ, tail_mean = tail_stats.mean, tail_sigma = tail_stats.sigma,
        # auxiliary pz parameters
        auxpz1_mean = auxpz1_stats.mean, auxpz1_sigma = auxpz1_stats.sigma, auxpz1_slope_sigma = auxpz1_stats.slope_residual_sigma,
        auxpz2_mean = auxpz2_stats.mean, auxpz2_sigma = auxpz2_stats.sigma, auxpz2_slope_sigma = auxpz2_stats.slope_residual_sigma,
        # timing parameters
        t0 = t0, t10 = t10, t50 = t50, t80 = t80, t90 = t90, t99 = t99, t50_pre = t50_pre,
        drift_time = drift_time, t50_current = t50_current, t_amax = t_amax, t_a_sg = t_a_sg, t_a_mwa = t_a_mwa,
        # energies (fixed and trap/cusp/zac) & extrema of trap/cusp/zac 
        e_10410 = e_10410_stats.max, e_535 = e_535_stats.max, e_313 = e_313_stats.max, e_101010 = e_101010_stats.max,
        t_10410 = e_10410_stats.tmax, t_535 = e_535_stats.tmax, t_313 = e_313_stats.tmax, t_101010 = e_101010_stats.tmax,
        e_trap = e_trap, e_cusp = e_cusp, e_zac = e_zac,
        e_trap_max = e_trap_extremestats.max, e_cusp_max = e_cusp_extremestats.max, e_zac_max = e_zac_extremestats.max,
        t_trap_max = e_trap_extremestats.tmax, t_cusp_max = e_cusp_extremestats.tmax, t_zac_max = e_zac_extremestats.tmax,
        # Qdrift & LQ
        qdrift = qdrift, lq80..., lq90..., lqamax...,
        # currents
        a_sg = a_sg, a_raw = a_raw, a_mwa = a_mwa, a_sg_fixed..., a_mwa_fixed...,
        # trigger
        trigger_position = triggers.intersect, trigger_multiplicity = triggers.n,
        # inverse parameters
        e_10410_inv = e_10410_max_inv, e_313_inv = e_313_max_inv,
        t0_inv = t0_inv,
    )
end
export dsp_icpc_compressed
