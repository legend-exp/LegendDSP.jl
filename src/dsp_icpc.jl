# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    dsp_icpc(data::Table, config::DSPConfig, τ::Quantity, pars_filter::PropDict; f_evaluate_qc::Union{Function, Missing}=missing)

DSP chain for HPGe (ICPC) detector waveforms. Returns a table with one row per waveform.

The baseline is subtracted, the waveforms are pole-zero corrected with the decay time `τ`, and timing, energy, current, charge-drift, quality and saturation parameters are determined. Use [`dsp_icpc_compressed`](@ref) for compressed data.

# Arguments
- `data`: table with the input columns below; encoded columns are decoded with `decode_data`
- `config`: [`DSPConfig`](@ref). `config.kwargs_pars` needs the entries `fc_bit_depth`, `t0_flt_pars`, `t0_mintot`, `tx_mintot`, `intrace_mintot`, `int_interpolation_order`, `int_interpolation_length`, `sig_interpolation_order` and `sig_interpolation_length`.
- `τ`: decay time of the preamplifier signal, used for the pole-zero correction
- `pars_filter`: optimized filter parameters as a `PropDict`: `trap`, `cusp` and `zac` with the entries `rt` and `ft`, and `sg` and `mwa` with the entry `wl`. Missing entries are taken from `config.default_flt_param`.

# Keywords
- `f_evaluate_qc`: function for [`get_qc_classifier`](@ref). With `missing` no classification is done and `qc_label` is -1.

# Input columns
`waveform`, `baseline`, `timestamp`, `eventnumber`, `daqenergy`

# Output columns
Amplitudes have the unit of the waveform samples. Times with a stated unit are converted; the others are in the unit of the waveform time axis.

Baseline and tail:
- `blmean`, `blsigma`, `blslope`, `bloffset`: statistics of the waveform in `config.bl_window`
- `tailmean`, `tailsigma`, `tailslope`, `tailoffset`: the same statistics of the pole-zero-corrected waveform in `config.tail_window`
- `tail_τ`, `tail_mean`, `tail_sigma`: decay time, and mean and standard deviation of the logarithm of the baseline-subtracted waveform in `config.tail_window`, from [`tailstats`](@ref)

Timing, from the pole-zero-corrected waveform, in µs:
- `t0`: start of the charge rise, from [`get_t0`](@ref)
- `t10`, `t50`, `t80`, `t90`, `t99`: first times at which the waveform reaches 10%, 50%, 80%, 90% and 99% of `e_101010`, from [`get_threshold`](@ref)
- `t50_current`: time at which the Savitzky-Golay derivative reaches half of its maximum
- `drift_time`: `t90 - t0`, in ns

Energies, from the pole-zero-corrected waveform:
- `e_max`, `e_min`: maximum and minimum of the baseline-subtracted waveform before the pole-zero correction
- `e_10410`, `e_535`, `e_313`, `e_101010`: maximum after a `TrapezoidalChargeFilter` with rise time and flat-top time of 10 and 4 µs, 5 and 3 µs, 3 and 1 µs, and 10 and 10 µs
- `e_trap`, `e_cusp`, `e_zac`: value of the trapezoidal, CUSP and ZAC filtered waveform at `t50 + rt + ft/2` (trapezoidal) and `t50 + flt_length/2` (CUSP and ZAC), determined by a signal estimator. `rt` and `ft` are the optimized rise and flat-top time of the respective filter in `pars_filter`, and `flt_length` is `config.flt_length_cusp` or `config.flt_length_zac`.
- `e_trap_max`, `e_cusp_max`, `e_zac_max`: maximum of the filtered waveform
- `t_trap_max`, `t_cusp_max`, `t_zac_max`: time of that maximum
- `e_10410_inv`, `e_313_inv`: maximum after the `e_10410` and `e_313` filters of the waveform multiplied by -1
- `t0_inv`: t0 of the waveform multiplied by -1, in µs, determined with the default `flt_pars` of [`get_t0`](@ref)

Current:
- `a_raw`: maximum in `config.current_window` of the simple derivative (`DerivativeFilter(1)`)
- `a_sg`, `t_a_sg`: maximum in `config.current_window` of the Savitzky-Golay derivative with window length `sg.wl` and polynomial degree `config.sg_flt_degree`, and its interpolated time
- `a_mwa`, `t_a_mwa`: maximum in `config.current_window` of the simple derivative after `MovingWindowMultiFilter(mwa.wl; passes = config.a_mwa_passes)`, and its interpolated time
- `a_sg_<ns>`, `a_mwa_<ns>`: the same maxima for each fixed window length in `config.a_sg_windows` and `config.a_mwa_windows`, with the window length in ns in the column name
- `t_amax`: interpolated time of the maximum of the simple derivative over the whole waveform

Charge drift:
- `qdrift`: Q-drift parameter from [`get_qdrift`](@ref), using `t0` and `config.qdrift_int_length`
- `lq80_<early><gap><late>`, `lq90_<early><gap><late>`, `lqamax_<early><gap><late>`: LQ parameters from [`get_lq`](@ref) starting at `t80`, `t90` and `t_amax`, one column for each `[early, gap, late]` of `config.lq_windows`. The window lengths in the name are in units of 0.1 µs with two digits each, so `lq80_250025` has an early window of 2.5 µs, no gap and a late window of 2.5 µs.

Quality and saturation:
- `qc_label`: label from `f_evaluate_qc`, evaluated on the baseline-subtracted waveform before the pole-zero correction
- `inTrace_intersect`, `inTrace_n`: time and number of in-trace pile-up triggers, from [`get_triggers`](@ref) on the Savitzky-Golay derivative with the threshold `config.inTraceCut_std_threshold` in standard deviations of its baseline
- `n_sat_low`, `n_sat_high`: number of samples of the unsubtracted waveform equal to 0 and to `2^fc_bit_depth - fc_bit_depth`, from [`saturation`](@ref)
- `n_sat_low_cons`, `n_sat_high_cons`: longest run of consecutive samples at those values

FADC values:
- `blfc`, `timestamp`, `eventID_fadc`, `e_fc`: baseline, timestamp, event number and energy from the FADC
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
    dsp_icpc_compressed(data::Table, config::DSPConfig, τ::Quantity, pars_filter::PropDict; f_evaluate_qc::Union{Function, Missing}=missing)

DSP chain for HPGe (ICPC) detector waveforms of compressed data. Returns a table with one row per waveform.

Each event has a presummed waveform, in which `presum_rate` consecutive samples are summed, and a windowed waveform with the original sampling around the signal rise. Energies, the tail and the baseline are determined from the presummed waveform, and timing, current and charge-drift parameters from the windowed waveform. The baseline mean of the presummed waveform is subtracted from it, and `blmean / presum_rate` from the windowed waveform. Both waveforms are then pole-zero corrected with `τ`.

# Arguments
- `data`: table with the input columns below; encoded columns are decoded with `decode_data`. `presum_rate` has to be the same for all events.
- `config`: [`DSPConfig`](@ref). `config.kwargs_pars` needs the entries `fc_bit_depth`, `t0_flt_pars`, `t0_mintot`, `tx_mintot`, `int_interpolation_order`, `int_interpolation_length`, `sig_interpolation_order` and `sig_interpolation_length`.
- `τ`: decay time of the preamplifier signal, used for the pole-zero correction
- `pars_filter`: optimized filter parameters as a `PropDict`, as for [`dsp_icpc`](@ref)

# Keywords
- `f_evaluate_qc`: function for [`get_qc_classifier_compressed`](@ref). With `missing` no classification is done and `qc_label` is -1.

# Input columns
`waveform_presummed`, `waveform_windowed`, `presum_rate`, `baseline`, `timestamp`, `eventnumber`, `daqenergy`, `t_sat_lo`, `t_sat_hi`, `deadtime`

# Output columns
Amplitudes have the unit of the waveform samples; those from the presummed waveform are sums over `presum_rate` samples. Times with a stated unit are converted; the others are in the unit of the waveform time axis.

FADC values and saturation:
- `blfc`, `timestamp`, `eventID_fadc`, `e_fc`, `deadtime`, `t_sat_lo`, `t_sat_hi`: copied from the input
- `n_sat_low`, `n_sat_high`, `n_sat_low_cons`, `n_sat_high_cons`: number of samples and longest run of consecutive samples of the unsubtracted presummed waveform equal to 0 and to `(2^fc_bit_depth - fc_bit_depth) * presum_rate`, from [`saturation`](@ref)

Baseline and tail, from the presummed waveform:
- `blmean`, `blsigma`, `blslope`, `bloffset`, `bl_slope_sigma`: statistics in `config.bl_window`, with `bl_slope_sigma` the standard deviation of the residuals of the slope fit
- `tailmean`, `tailsigma`, `tailslope`, `tailoffset`: the same statistics of the pole-zero-corrected waveform in `config.tail_window`
- `tail_τ`, `tail_mean`, `tail_sigma`: decay time, and mean and standard deviation of the logarithm of the baseline-subtracted waveform in `config.tail_window`, from [`tailstats`](@ref)
- `auxpz1_mean`, `auxpz1_sigma`, `auxpz1_slope_sigma`, `auxpz2_mean`, `auxpz2_sigma`, `auxpz2_slope_sigma`: mean, standard deviation and slope residual of the pole-zero-corrected waveform in `config.auxpz1_window` and `config.auxpz2_window`

Extrema:
- `e_max`, `t_max`: maximum and its time from the baseline-subtracted presummed waveform
- `e_min`, `t_min`: minimum and its time from the baseline-subtracted windowed waveform
- `e_max_pre`, `e_min_pre`, `t_max_pre`, `t_min_pre`: maximum, minimum and their times from the baseline-subtracted presummed waveform

Timing, in µs, from the pole-zero-corrected waveforms:
- `t0`: start of the charge rise of the windowed waveform, from [`get_t0`](@ref)
- `t10`, `t50`, `t80`, `t90`, `t99`: first times at which the windowed waveform reaches 10%, 50%, 80%, 90% and 99% of `e_101010 / presum_rate`
- `t50_pre`: first time at which the presummed waveform reaches 50% of `e_101010`
- `drift_time`: `t90 - t0`, in ns
- `t50_current`: time at which the Savitzky-Golay derivative of the presummed waveform reaches half of its maximum. The filter window length is `sg.wl * presum_rate / 2`.
- `t_amax`: interpolated time of the maximum of the simple derivative of the windowed waveform

Energies, from the presummed waveform:
- `e_10410`, `e_535`, `e_313`, `e_101010`: maximum after a `TrapezoidalChargeFilter` with rise time and flat-top time of 10 and 4 µs, 5 and 3 µs, 3 and 1 µs, and 10 and 10 µs
- `t_10410`, `t_535`, `t_313`, `t_101010`: times of those maxima
- `e_trap`, `e_cusp`, `e_zac`: value of the trapezoidal, CUSP and ZAC filtered waveform at `t50_pre + rt + ft/2` (trapezoidal) and `t50_pre + flt_length/2` (CUSP and ZAC), determined by a signal estimator, with the parameters as for [`dsp_icpc`](@ref)
- `e_trap_max`, `e_cusp_max`, `e_zac_max`, `t_trap_max`, `t_cusp_max`, `t_zac_max`: maximum of the filtered waveform and its time
- `e_10410_inv`, `e_313_inv`: maximum after the `e_10410` and `e_313` filters of the presummed waveform multiplied by -1
- `t0_inv`: t0 of the windowed waveform multiplied by -1, in µs, determined with the default `flt_pars` of [`get_t0`](@ref)

Current, from the windowed waveform:
- `a_raw`, `a_sg`, `t_a_sg`, `a_mwa`, `t_a_mwa`, `a_sg_<ns>`, `a_mwa_<ns>`: as for [`dsp_icpc`](@ref)

Charge drift, from the windowed waveform:
- `qdrift`: Q-drift parameter from [`get_qdrift`](@ref), using `t0` and `config.qdrift_int_length`
- `lq80_<early><gap><late>`, `lq90_<early><gap><late>`, `lqamax_<early><gap><late>`: LQ parameters as for [`dsp_icpc`](@ref)

Quality:
- `qc_label`: label from `f_evaluate_qc`, evaluated on the baseline-subtracted presummed waveform before the pole-zero correction
- `trigger_position`, `trigger_multiplicity`: time and number of in-trace triggers, from [`get_triggers`](@ref) on the presummed waveform filtered with `TrapezoidalChargeFilter(2u"µs", 1u"µs")`, with the threshold `config.inTraceCut_std_threshold` in standard deviations of its baseline and a minimum time over threshold of 2 µs
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
