
"""
    dsp_trap_rt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}) where T<:Real

Get ENC noise grid values for given trap grid rise times.

Returns:
    - `enc_trap_grid`: Array ENC noise values for the given trap rise time grid
"""
function dsp_trap_rt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, ft::Quantity{T}=4.0u"µs") where T<:Real
    # get config parameters
    bl_mean_min, bl_mean_max    = config.bl_mean
    e_grid_rt_trap              = config.e_grid_rt_trap
    t0_threshold                = config.t0_threshold
    enc_pickoff_trap            = config.enc_pickoff_trap

    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs, bl_mean_min, bl_mean_max)

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)

    # deconvolute waveform
    deconv_flt = InvCRFilter(τ)
    wvfs_pz = deconv_flt.(wvfs)

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

    # truncate waveform for ENC filtering
    uflt_trunc_enc = TruncateFilter(0u"µs"..40u"µs")
    wvfs_pz_trunc = uflt_trunc_enc.(wvfs_pz)

    # get energy grid for efficient optimization
    enc_trap_grid = zeros(Float64, length(e_grid_rt_trap), length(wvfs_pz))
    for (r, rt) in enumerate(e_grid_rt_trap)
        if rt < ft
            continue
        end
        uflt_rtft      = TrapezoidalChargeFilter(rt, ft)
        
        wvfs_flt_rtft  = uflt_rtft.(wvfs_pz_trunc)

        enc_rtft       = SignalEstimator(PolynomialDNI(4, 80u"ns")).(wvfs_flt_rtft, enc_pickoff_trap)

        enc_trap_grid[r, :]   = enc_rtft
    end

    return enc_trap_grid
end
export dsp_trap_rt_optimization


"""
    dsp_trap_ft_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, rt::Quantity{T}) where T<:Real

Get energy grid values for given trap grid rise times while varying the flat-top time.

Returns:
    - `e_grid`: Array energy values for the given trap rise time grid at a given flat-top time grid.
"""
function dsp_trap_ft_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, rt::Quantity{T}) where T<:Real
    # get config parameters
    bl_mean_min, bl_mean_max    = config.bl_mean
    e_grid_ft_trap              = config.e_grid_ft_trap
    t0_threshold                = config.t0_threshold
    enc_pickoff_trap            = config.enc_pickoff_trap

    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs, bl_mean_min, bl_mean_max)

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)

    # deconvolute waveform
    deconv_flt = InvCRFilter(τ)
    wvfs_pz = deconv_flt.(wvfs)

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

    # get energy grid for efficient optimization
    e_grid   = Array{Union{Missing, Float32}}(missing, length(e_grid_ft_trap), length(wvfs_pz))
    for (f, ft) in enumerate(e_grid_ft_trap)
        if rt < ft
            continue
        end
        uflt_rtft      = TrapezoidalChargeFilter(rt, ft)
        
        wvfs_flt_rtft  = uflt_rtft.(wvfs_pz)

        e_rtft         = SignalEstimator(PolynomialDNI(4, 80u"ns")).(wvfs_flt_rtft, t0 .+ (rt + ft/2))

        e_grid[f, :]     = e_rtft
    end

    return e_grid
end
export dsp_trap_ft_optimization


""" 
    dsp_sg_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, pars_filter::PropDict) where T<:Real

Optimize the Savitzky-Golay filter parameters for a given waveform set.

Returns:
    - `aoe`: Array of efficiency values for the given Savitzky-Golay filter parameters
    - `e`: Array of energy values for the given Savitzky-Golay filter parameters
    - `blmean`: Baseline mean value
    - `blslope`: Baseline slope value
"""
function dsp_sg_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, pars_filter::PropDict) where T<:Real
    # get config parameters
    bl_mean_min, bl_mean_max    = config.bl_mean
    a_grid_wl_sg                = config.a_grid_wl_sg
    t0_threshold                = config.t0_threshold
    enc_pickoff_trap            = config.enc_pickoff_trap

    # get optimal filter parameters
    rt = pars_filter.trap_rt.val*u"µs"
    ft = pars_filter.trap_ft.val*u"µs"

    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs, bl_mean_min, bl_mean_max)

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)

    # deconvolute waveform
    deconv_flt = InvCRFilter(τ)
    wvfs_pz = deconv_flt.(wvfs)

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

    # get energy for filter parameters
    uflt_rtft      = TrapezoidalChargeFilter(rt, ft)
    wvfs_flt_rtft  = uflt_rtft.(wvfs_pz)
    e_rtft = SignalEstimator(PolynomialDNI(4, 80u"ns")).(wvfs_flt_rtft, t0 .+ (rt + ft/2))

    # extract current with filter length in grid with second order polynominal and first derivative
    aoe_grid   = Array{Union{Missing, Float32}}(missing, length(a_grid_wl_sg), length(wvfs_pz))
    for (w, wl) in enumerate(a_grid_wl_sg)
        sgflt_deriv = SavitzkyGolayFilter(wl, 2, 1)
        wvfs_sgflt_deriv = sgflt_deriv.(wvfs_pz)
        current_max = get_wvf_maximum.(wvfs_sgflt_deriv, 20u"µs", 100u"µs")

        aoe_grid[w, :]     = ustrip.(current_max) ./ e_rtft
    end
    return (aoe = aoe_grid, e = e_rtft, blmean = bl_stats.mean, blslope = bl_stats.slope)
end
export dsp_sg_optimization
