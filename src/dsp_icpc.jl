

function dsp_icpc(data::Q, config::DSPConfig, τ::Quantity{T}, pars_filter::PropDict) where {Q <: Table, T<:Real}
    # get config parameters
    bl_mean_min, bl_mean_max    = config.bl_mean
    t0_threshold                = config.t0_threshold
    pz_fit_min, pz_fit_max      = config.pz_fit
    inTraceCut_std_threshold    = config.inTraceCut_std_threshold

    # get optimal filter parameters
    trap_rt = pars_filter.trap_rt.val*u"µs"
    trap_ft = pars_filter.trap_ft.val*u"µs"
    sg_wl   = pars_filter.sg_wl.val*u"ns"


    # get waveform data 
    wvfs = data.waveform
    blfc = data.baseline
    ts   = data.timestamp
    evID = data.eventnumber
    efc  = data.daqenergy

    # get number of samples the waveform is saturated at low and high of FADC range
    bit_depth = 16 # of FlashCam FADC
    sat_low, sat_high = 0, 2^bit_depth - bit_depth
    sat_stats = saturation.(wvfs, sat_low, sat_high)

    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs, bl_mean_min, bl_mean_max)

    # pretrace difference 
    pretrace_diff = flatview(wvfs.signal)[1, :] - bl_stats.mean

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)

    # extract decay times
    tail_stats = tailstats.(wvfs, pz_fit_min, pz_fit_max)

    # deconvolute waveform
    deconv_flt = InvCRFilter(τ)
    wvfs_pz = deconv_flt.(wvfs)

    # get tail mean, std and slope
    pz_stats = signalstats.(wvfs_pz, pz_fit_min, pz_fit_max)

    # t0 determination
    # filter with fast asymetric trapezoidal filter and truncate waveform
    uflt_asy_t0 = TrapezoidalChargeFilter(40u"ns", 100u"ns", 2000u"ns")
    uflt_trunc_t0 = TruncateFilter(0u"µs"..60u"µs")

    # eventuell zwei schritte!!!
    wvfs_flt_asy_t0 = uflt_asy_t0.(uflt_trunc_t0.(wvfs_pz))

    # get intersect at t0 threshold (fixed as in MJD analysis)
    flt_intersec_t0 = Intersect(mintot = 600u"ns")

    # get t0 for every waveform as pick-off at fixed threshold
    t0 = uconvert.(u"µs", flt_intersec_t0.(wvfs_flt_asy_t0, t0_threshold).x)

    # get risetimes and drift times by intersection
    flt_intersec_90RT = Intersect(mintot = 100u"ns")
    flt_intersec_99RT = Intersect(mintot = 20u"ns")
    flt_intersec_lowRT = Intersect(mintot = 600u"ns")
    
    wvf_max = maximum.(wvfs.signal)

    rt1090     = uconvert.(u"ns", flt_intersec_90RT.(wvfs_pz, wvf_max .* 0.9).x - flt_intersec_lowRT.(wvfs_pz, wvf_max .* 0.1).x)
    rt1099     = uconvert.(u"ns", flt_intersec_99RT.(wvfs_pz, wvf_max .* 0.99).x - flt_intersec_lowRT.(wvfs_pz, wvf_max .* 0.1).x)
    rt9099     = uconvert.(u"ns", flt_intersec_99RT.(wvfs_pz, wvf_max .* 0.99).x - flt_intersec_90RT.(wvfs_pz, wvf_max .* 0.90).x)
    drift_time = uconvert.(u"ns", flt_intersec_90RT.(wvfs_pz, wvf_max .* 0.90).x - t0)

    # get Q-drift parameter
    int_flt = IntegratorFilter(1)
    wvfs_flt_int = int_flt.(wvfs_pz)

    area1 = SignalEstimator(PolynomialDNI(3, 100u"ns")).(wvfs_flt_int, t0 .+ 2.5u"µs") .- SignalEstimator(PolynomialDNI(3, 100u"ns")).(wvfs_flt_int, t0)
    area2 = SignalEstimator(PolynomialDNI(3, 100u"ns")).(wvfs_flt_int, t0 .+ 5u"µs")   .- SignalEstimator(PolynomialDNI(3, 100u"ns")).(wvfs_flt_int, t0 .+ 2.5u"µs")
    qdrift = area2 .- area1

    # extract energy and ENC noise param from maximum of filtered wvfs
    uflt_10410 = TrapezoidalChargeFilter(10u"µs", 4u"µs")

    wvfs_flt_10410 = uflt_10410.(wvfs_pz)
    e_10410        = SignalEstimator(PolynomialDNI(3, 100u"ns")).(wvfs_flt_10410, t0 .+ 12u"µs")

    uflt_zac10410 = ZACChargeFilter(10u"µs", 4u"µs", 30u"µs")

    wvfs_flt_zac10410 = uflt_zac10410.(wvfs_pz)
    e_zac_10410       = SignalEstimator(PolynomialDNI(3, 100u"ns")).(wvfs_flt_zac10410, t0 .+ 12u"µs")


    # get energy of optimized rise and flat-top time
    uflt_rtft = TrapezoidalChargeFilter(trap_rt, trap_ft)

    wvfs_flt_rtft  = uflt_rtft.(wvfs_pz)
    e_rtft         = SignalEstimator(PolynomialDNI(3, 100u"ns")).(wvfs_flt_rtft, t0 .+ (trap_rt + trap_ft/2))


    # extract current with filter length of 180ns with second order polynominal and first derivative
    sgflt_deriv = SavitzkyGolayFilter(sg_wl, 2, 1)
    wvfs_sgflt_deriv = sgflt_deriv.(wvfs_pz)
    current_max = get_wvf_maximum.(wvfs_sgflt_deriv, 20u"µs", 100u"µs")

    # in-trace pile-up rejector
    flt_intersec_inTrace = Intersect(mintot = 100u"ns")
    deriv_stats = signalstats.(wvfs_sgflt_deriv, bl_mean_min + wvfs_sgflt_deriv.time[1][1], bl_mean_max)
    # threshold over std
    inTraceCut = inTraceCut_std_threshold .* deriv_stats.sigma

    # get position and multiplicity of in-trace pile-up
    inTrace_pileUp      = flt_intersec_inTrace.(reverse_waveform.(wvfs_sgflt_deriv), inTraceCut)
    inTrace_intersect   = wvfs_pz.time[1][end] .- inTrace_pileUp.x
    inTrace_n           = inTrace_pileUp.multiplicity



    # output Table 
    return TypedTables.Table(blmean = bl_stats.mean, blsigma = bl_stats.sigma, blslope = bl_stats.slope, bloffset = bl_stats.offset, 
    tailmean = pz_stats.mean, tailsigma = pz_stats.sigma, tailslope = pz_stats.slope, tailoffset = pz_stats.offset,
    t0 = t0, tail_τ = tail_stats.τ, tail_mean = tail_stats.mean, tail_sigma = tail_stats.sigma,
    e_10410 = e_10410, e_zac_10410 = e_zac_10410,
    e_trap = e_rtft,
    qdrift = qdrift,
    a = current_max,
    blfc = blfc, timestamp = ts, eventID_fadc = evID, e_fc = efc,
    pretrace_diff = pretrace_diff, 
    rt1090 = rt1090, rt1099 = rt1099, rt9099 = rt9099, drift_time = drift_time,
    inTrace_intersect = inTrace_intersect, inTrace_n = inTrace_n,
    n_sat_low = sat_stats.low, n_sat_high = sat_stats.high, n_sat_low_cons = sat_stats.max_cons_low, n_sat_high_cons = sat_stats.max_cons_high
    )
end
export dsp_icpc