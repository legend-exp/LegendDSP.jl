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
- `t10`: timepoint of 10% of waveform maximum
- `t50`: timepoint of 50% of waveform maximum
- `t80`: timepoint of 80% of waveform maximum
- `t90`: timepoint of 90% of waveform maximum
- `t99`: timepoint of 99% of waveform maximum
- `t50_current`: timepoint of current rise to 50% of maximum
- `tail_τ`: tail decay time
- `tail_mean`: tail mean before PZ correction
- `tail_sigma`: tail sigma before PZ correction
- `e_max`: maximum of waveform
- `e_min`: minimum of waveform
- `e_10410`: energy of waveform with trapezoidal filter of 10µs rise time with 4µs flat-top
- `e_313`: energy of waveform with trapezoidal filter of 3µs rise time with 1µs flat-top
- `e_10410_inv`: maximum of inverted waveform with trapezoidal filter of 10µs rise time with 4µs flat-top
- `e_313_inv`: maximum of inverted waveform with trapezoidal filter of 3µs rise time with 1µs flat-top
- `t0_inv`: start time of inverted waveform drift
- `e_trap`: energy of waveform with trapezoidal filter of optimized rise and flat-top time
- `e_cusp`: energy of waveform with CUSP filter of optimized rise and flat-top time
- `e_zac`: energy of waveform with ZAC filter of optimized rise and flat-top time
- `qdrift`: Q-drift parameter
- `lq`: LQ parameter
- `a`: current maximum with optimal Savitzky-Golay filter length parameter
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
function dsp_icpc(data::Q, config::DSPConfig, τ::Quantity{T}, pars_filter::PropDict) where {Q <: Table, T<:Real}
    # get config parameters
    bl_mean_min, bl_mean_max    = config.bl_mean
    t0_threshold                = config.t0_threshold
    pz_fit_min, pz_fit_max      = config.pz_fit
    inTraceCut_std_threshold    = config.inTraceCut_std_threshold
    
    # get optimal filter parameters
    trap_rt = pars_filter.trap.rt.val*u"µs"
    trap_ft = pars_filter.trap.ft.val*u"µs"
    cusp_rt = pars_filter.cusp.rt.val*u"µs"
    cusp_ft = pars_filter.cusp.ft.val*u"µs"
    zac_rt  = pars_filter.zac.rt.val*u"µs"
    zac_ft  = pars_filter.zac.ft.val*u"µs"
    sg_wl   = pars_filter.sg.wl.val*u"ns"

    # get waveform data 
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
    bit_depth = 16 # of FlashCam FADC
    sat_low, sat_high = 0, 2^bit_depth - bit_depth
    sat_stats = saturation.(wvfs, sat_low, sat_high)

    # set tau for CUSP filter to very high number to switch of CR filter
    τ_cusp = 10000000.0u"µs"
    τ_zac = 10000000.0u"µs"

    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs, bl_mean_min, bl_mean_max)

    # pretrace difference 
    pretrace_diff = flatview(wvfs.signal)[1, :] - bl_stats.mean

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)

    # get wvf maximum
    wvf_max = maximum.(wvfs.signal)
    wvf_min = minimum.(wvfs.signal)

    # extract decay times
    tail_stats = tailstats.(wvfs, pz_fit_min, pz_fit_max)

    # deconvolute waveform 
    # --> wvfs = wvfs_pz
    deconv_flt = InvCRFilter(τ)
    wvfs = deconv_flt.(wvfs)

    # get tail mean, std and slope
    pz_stats = signalstats.(wvfs, pz_fit_min, pz_fit_max)

    # t0 determination
    t0 = get_t0(wvfs, t0_threshold)

    # get threshold points in rise
    t10 = get_threshold(wvfs, wvf_max .* 0.1)
    t50 = get_threshold(wvfs, wvf_max .* 0.5)
    t80 = get_threshold(wvfs, wvf_max .* 0.8)
    t90 = get_threshold(wvfs, wvf_max .* 0.9)
    t99 = get_threshold(wvfs, wvf_max .* 0.99)
    
    drift_time = uconvert.(u"ns", t90 - t0)

    # get Q-drift parameter
    qdrift = get_qdrift(wvfs, t0, (2.5:5.0)u"µs")

    # get LQ parameter
    lq  = get_qdrift(wvfs, t80, (2.5:5.0)u"µs")

    # robust energy reconstruction with long, middle and short rise and flat-top times
    uflt_10410 = TrapezoidalChargeFilter(10u"µs", 4u"µs")
    e_10410  = maximum.((uflt_10410.(wvfs)).signal)

    uflt_535 = TrapezoidalChargeFilter(5u"µs", 3u"µs")
    e_535  = maximum.((uflt_535.(wvfs)).signal)

    uflt_313 = TrapezoidalChargeFilter(3u"µs", 1u"µs")
    e_313  = maximum.((uflt_313.(wvfs)).signal)

    # signal estimator for precise energy reconstruction
    signal_estimator = SignalEstimator(PolynomialDNI(3, 100u"ns"))

    # get trap energy of optimized rise and flat-top time
    uflt_trap_rtft = TrapezoidalChargeFilter(trap_rt, trap_ft)

    e_trap = signal_estimator.(uflt_trap_rtft.(wvfs), t50 .+ (trap_rt + trap_ft/2))

    # get cusp energy of optimized rise and flat-top time
    uflt_cusp_rtft = CUSPChargeFilter(cusp_rt, cusp_ft, τ_cusp, flt_length_cusp, cusp_scale)

    e_cusp = signal_estimator.(uflt_cusp_rtft.(wvfs), t50 .+ (flt_length_cusp /2))

    # get zac energy of optimized rise and flat-top time
    uflt_zac_rtft = ZACChargeFilter(zac_rt, zac_ft, τ_zac, flt_length_zac, zac_scale)

    e_zac = signal_estimator.(uflt_zac_rtft.(wvfs), t50 .+ (flt_length_zac /2))

    # extract current with optimal SG filter length with second order polynominal and first derivative
    wvfs_sgflt_deriv = SavitzkyGolayFilter(sg_wl, 2, 1).(wvfs)

    current_window = (20.0:100.0)u"µs"
    current_max = get_wvf_maximum.(wvfs_sgflt_deriv, first(current_window), last(current_window))

    # get in-trace pile-up
    inTrace_pileUp = get_intracePileUp(wvfs_sgflt_deriv, inTraceCut_std_threshold, bl_mean_min, bl_mean_max)
    
    # get position of current rise
    t50_current = get_threshold(wvfs_sgflt_deriv, maximum.(wvfs_sgflt_deriv.signal) .* 0.5; mintot=300u"ns")

    # invert waveform for DC tagging
    # wvfs --> wvfs_pz_inv
    wvfs = multiply_waveform.(wvfs, -1.0)

    # get inverted waveform maximum for long and short filter
    e_10410_max_inv  = maximum.(uflt_10410.(wvfs).signal)

    e_313_max_inv  = maximum.(uflt_313.(wvfs).signal)

    # t0 determination
    t0_inv = get_t0(wvfs, t0_threshold)

    # output Table 
    return TypedTables.Table(blmean = bl_stats.mean, blsigma = bl_stats.sigma, blslope = bl_stats.slope, bloffset = bl_stats.offset, 
    tailmean = pz_stats.mean, tailsigma = pz_stats.sigma, tailslope = pz_stats.slope, tailoffset = pz_stats.offset,
    t0 = t0, t10 = t10, t50 = t50, t80 = t80, t90 = t90, t99 = t99,
    t50_current = t50_current, 
    drift_time = drift_time,
    tail_τ = tail_stats.τ, tail_mean = tail_stats.mean, tail_sigma = tail_stats.sigma,
    e_max = wvf_max, e_min = wvf_min,
    e_10410 = e_10410, e_535 = e_535, e_313 = e_313,
    e_10410_inv = e_10410_max_inv, e_313_inv = e_313_max_inv,
    t0_inv = t0_inv,
    e_trap = e_trap, e_cusp = e_cusp, e_zac = e_zac, 
    qdrift = qdrift, lq = lq,
    a = current_max,
    blfc = blfc, timestamp = ts, eventID_fadc = evID, e_fc = efc,
    pretrace_diff = pretrace_diff, 
    inTrace_intersect = inTrace_pileUp.intersect, inTrace_n = inTrace_pileUp.n,
    n_sat_low = sat_stats.low, n_sat_high = sat_stats.high, n_sat_low_cons = sat_stats.max_cons_low, n_sat_high_cons = sat_stats.max_cons_high
    )
end
export dsp_icpc