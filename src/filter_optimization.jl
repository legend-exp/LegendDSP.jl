
"""
    dsp_trap_rt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}) where T<:Real

Get ENC noise grid values for given trap grid rise times.

Returns:
    - `enc_trap_grid`: Array ENC noise values for the given trap rise time grid
"""
function dsp_trap_rt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T},; ft::Quantity{T}=4.0u"µs", flt_length::Quantity{T}=20.0u"µs") where T<:Real
    # get config parameters
    bl_mean_min, bl_mean_max    = config.bl_mean
    e_grid_rt_trap              = config.e_grid_rt_trap
    enc_pickoff_trap            = config.enc_pickoff_trap

    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs, bl_mean_min, bl_mean_max)

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)

    # deconvolute waveform
    deconv_flt = InvCRFilter(τ)
    wvfs_pz = deconv_flt.(wvfs)

    # get energy grid for efficient optimization
    enc_trap_grid = zeros(Float64, length(e_grid_rt_trap), length(wvfs_pz))
    for (r, rt) in enumerate(e_grid_rt_trap)
        if rt < ft
            continue
        end
        uflt_rtft      = TrapezoidalChargeFilter(rt, ft)
        
        wvfs_flt_rtft  = uflt_rtft.(wvfs_pz)

        enc_rtft       = SignalEstimator(PolynomialDNI(4, 80u"ns")).(wvfs_flt_rtft, enc_pickoff_trap)

        enc_trap_grid[r, :]   = enc_rtft
    end

    return enc_trap_grid
end
export dsp_trap_rt_optimization


"""
    dsp_cusp_rt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T},; ft::Quantity{T}=4.0u"µs", length::Quantity{T}=10.0u"µs") where T<:Real

Get ENC noise grid values for given CUSP grid rise times.

Returns:
    - `enc_cusp_grid`: Array ENC noise values for the given CUSP rise time grid
"""
function dsp_cusp_rt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T},; ft::Quantity{T}=4.0u"µs") where T<:Real
    # get config parameters
    bl_mean_min, bl_mean_max    = config.bl_mean
    e_grid_rt_cusp              = config.e_grid_rt_cusp
    enc_pickoff_cusp            = config.enc_pickoff_cusp
    flt_length_cusp             = config.flt_length_cusp
    cusp_scale                  = ustrip(NoUnits, flt_length_cusp/step(wvfs[1].time))

    # set tau for CUSP filter to very high number to switch of CR filter
    τ_cusp = 10000000.0u"µs"

    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs, bl_mean_min, bl_mean_max)

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)

    # # truncate waveform for ENC filtering
    # uflt_trunc_enc = TruncateFilter(0u"µs"..40u"µs")
    # wvfs_trunc = uflt_trunc_enc.(wvfs)

    # deconvolute waveform
    deconv_flt = InvCRFilter(τ)
    wvfs_pz = deconv_flt.(wvfs)

    # get energy grid for efficient optimization
    enc_cusp_grid = zeros(Float64, length(collect(e_grid_rt_cusp)), length(wvfs))
    for (r, rt) in enumerate(e_grid_rt_cusp)
        if rt < ft
            continue
        end
        uflt_rtft      = CUSPChargeFilter(rt, ft, τ_cusp, flt_length_cusp, cusp_scale)
        
        wvfs_flt_rtft  = uflt_rtft.(wvfs_pz)

        enc_rtft       = SignalEstimator(PolynomialDNI(4, 80u"ns")).(wvfs_flt_rtft, enc_pickoff_cusp)

        enc_cusp_grid[r, :]   = enc_rtft
    end

    return enc_cusp_grid
end
export dsp_cusp_rt_optimization


"""
    dsp_zac_rt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T},; ft::Quantity{T}=4.0u"µs", length::Quantity{T}=10.0u"µs") where T<:Real

Get ENC noise grid values for given ZAC grid rise times.

Returns:
    - `enc_zac_grid`: Array ENC noise values for the given ZAC rise time grid
"""
function dsp_zac_rt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T},; ft::Quantity{T}=4.0u"µs") where T<:Real
    # get config parameters
    bl_mean_min, bl_mean_max    = config.bl_mean
    e_grid_rt_zac               = config.e_grid_rt_zac
    enc_pickoff_zac             = config.enc_pickoff_zac
    flt_length_zac              = config.flt_length_zac
    zac_scale                   = ustrip(NoUnits, flt_length_zac/step(wvfs[1].time))

    # set tau for ZAC filter to very high number to switch of CR filter
    τ_zac = 10000000.0u"µs"

    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs, bl_mean_min, bl_mean_max)

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)

    # truncate waveform for ENC filtering
    # uflt_trunc_enc = TruncateFilter(0u"µs"..40u"µs")
    # wvfs_trunc = uflt_trunc_enc.(wvfs)

    # deconvolute waveform
    deconv_flt = InvCRFilter(τ)
    wvfs_pz = deconv_flt.(wvfs)

    # get energy grid for efficient optimization
    enc_zac_grid = zeros(Float64, length(e_grid_rt_zac), length(wvfs))
    for (r, rt) in enumerate(e_grid_rt_zac)
        if rt < ft
            continue
        end
        uflt_rtft      = ZACChargeFilter(rt, ft, τ_zac, flt_length_zac, zac_scale)
        
        wvfs_flt_rtft  = uflt_rtft.(wvfs_pz)

        enc_rtft       = SignalEstimator(PolynomialDNI(4, 80u"ns")).(wvfs_flt_rtft, enc_pickoff_zac)

        enc_zac_grid[r, :]   = enc_rtft
    end

    return enc_zac_grid
end
export dsp_zac_rt_optimization


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

    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs, bl_mean_min, bl_mean_max)

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)

    # deconvolute waveform
    deconv_flt = InvCRFilter(τ)
    wvfs_pz = deconv_flt.(wvfs)

    # t50 determination
    wvf_max = maximum.(wvfs.signal)
    t50 = get_t50(wvfs_pz, wvf_max)

    # get energy grid for efficient optimization
    e_grid   = Array{Union{Missing, Float32}}(missing, length(e_grid_ft_trap), length(wvfs_pz))
    for (f, ft) in enumerate(e_grid_ft_trap)
        if rt < ft
            continue
        end
        uflt_rtft      = TrapezoidalChargeFilter(rt, ft)
        
        wvfs_flt_rtft  = uflt_rtft.(wvfs_pz)

        e_rtft         = SignalEstimator(PolynomialDNI(4, 80u"ns")).(wvfs_flt_rtft, t50 .+ (rt + ft/2))

        e_grid[f, :]     = e_rtft
    end

    return e_grid
end
export dsp_trap_ft_optimization


"""
    dsp_cusp_ft_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, rt::Quantity{T},; flt_length::Quantity{T}=10.0u"µs") where T<:Real

Get energy grid values for given CUSP grid rise times while varying the flat-top time.

Returns:
    - `e_grid`: Array energy values for the given CUSP rise time grid at a given flat-top time grid.
"""
function dsp_cusp_ft_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, rt::Quantity{T}) where T<:Real
    # get config parameters
    bl_mean_min, bl_mean_max    = config.bl_mean
    e_grid_ft_cusp              = config.e_grid_ft_cusp
    flt_length_cusp             = config.flt_length_cusp
    cusp_scale                  = ustrip(NoUnits, flt_length_cusp/step(wvfs[1].time))

    # set tau for CUSP filter to very high number to switch of CR filter
    τ_cusp = 10000000.0u"µs"

    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs, bl_mean_min, bl_mean_max)

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)

    # deconvolute waveform
    deconv_flt = InvCRFilter(τ)
    wvfs_pz = deconv_flt.(wvfs)

    # t50 determination
    wvf_max = maximum.(wvfs.signal)
    t50 = get_t50(wvfs_pz, wvf_max)

    # get energy grid for efficient optimization
    e_grid   = Array{Union{Missing, Float32}}(missing, length(e_grid_ft_cusp), length(wvfs))
    for (f, ft) in enumerate(e_grid_ft_cusp)
        if rt < ft
            continue
        end
        uflt_rtft      = CUSPChargeFilter(rt, ft, τ_cusp, flt_length_cusp, cusp_scale)
        
        wvfs_flt_rtft  = uflt_rtft.(wvfs_pz)

        e_rtft         = SignalEstimator(PolynomialDNI(4, 80u"ns")).(wvfs_flt_rtft, t50 .+ flt_length_cusp/2)

        e_grid[f, :]     = e_rtft
    end

    return e_grid
end
export dsp_cusp_ft_optimization


"""
    dsp_zac_ft_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, rt::Quantity{T},; flt_length::Quantity{T}=10.0u"µs") where T<:Real

Get energy grid values for given ZAC grid rise times while varying the flat-top time.

Returns:
    - `e_grid`: Array energy values for the given ZAC rise time grid at a given flat-top time grid.
"""
function dsp_zac_ft_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, rt::Quantity{T}) where T<:Real
    # get config parameters
    bl_mean_min, bl_mean_max    = config.bl_mean
    e_grid_ft_zac               = config.e_grid_ft_zac
    flt_length_zac              = config.flt_length_zac
    zac_scale                   = ustrip(NoUnits, flt_length_zac/step(wvfs[1].time))

    # set tau for ZAC filter to very high number to switch of CR filter
    τ_zac = 10000000.0u"µs"

    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs, bl_mean_min, bl_mean_max)

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)

    # deconvolute waveform
    deconv_flt = InvCRFilter(τ)
    wvfs_pz = deconv_flt.(wvfs)

    # t50 determination
    wvf_max = maximum.(wvfs.signal)
    t50 = get_t50(wvfs_pz, wvf_max)

    # get energy grid for efficient optimization
    e_grid   = Array{Union{Missing, Float32}}(missing, length(e_grid_ft_zac), length(wvfs))
    for (f, ft) in enumerate(e_grid_ft_zac)
        if rt < ft
            continue
        end
        uflt_rtft      = ZACChargeFilter(rt, ft, τ_zac, flt_length_zac, zac_scale)
        
        wvfs_flt_rtft  = uflt_rtft.(wvfs_pz)

        e_rtft         = SignalEstimator(PolynomialDNI(4, 80u"ns")).(wvfs_flt_rtft, t50 .+ flt_length_zac/2)

        e_grid[f, :]     = e_rtft
    end

    return e_grid
end
export dsp_zac_ft_optimization


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

    # get optimal filter parameters
    rt = pars_filter.trap.rt.val*u"µs"
    ft = pars_filter.trap.ft.val*u"µs"

    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs, bl_mean_min, bl_mean_max)

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)

    # deconvolute waveform
    deconv_flt = InvCRFilter(τ)
    wvfs_pz = deconv_flt.(wvfs)
    
    # t50 determination
    wvf_max = maximum.(wvfs.signal)
    t50 = get_t50(wvfs_pz, wvf_max)

    # get energy for filter parameters
    uflt_rtft      = TrapezoidalChargeFilter(rt, ft)
    wvfs_flt_rtft  = uflt_rtft.(wvfs_pz)
    e_rtft = SignalEstimator(PolynomialDNI(4, 80u"ns")).(wvfs_flt_rtft, t50 .+ (rt + ft/2))

    # extract current with filter length in grid with second order polynominal and first derivative
    aoe_grid   = ones(Float64, length(a_grid_wl_sg), length(wvfs_pz))
    for (w, wl) in enumerate(a_grid_wl_sg)
        sgflt_deriv = SavitzkyGolayFilter(wl, 2, 1)
        wvfs_sgflt_deriv = sgflt_deriv.(wvfs_pz)
        current_max = get_wvf_maximum.(wvfs_sgflt_deriv, 20u"µs", 100u"µs")

        aoe_grid[w, :]     = ustrip.(current_max) ./ e_rtft
    end
    return (aoe = aoe_grid, e = e_rtft, blmean = bl_stats.mean, blslope = bl_stats.slope, t50 = t50)
end
export dsp_sg_optimization
