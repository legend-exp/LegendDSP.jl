# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    dsp_trap_rt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T},; ft::Quantity{T}=4.0u"µs") where T<:Real

Get ENC noise grid values for given trap grid rise times.

# Returns
    - `enc_trap_grid`: Array ENC noise values for the given trap rise time grid
"""
function dsp_trap_rt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T},; ft::Quantity{T}=2.0u"µs") where T<:Real
    # get config parameters
    bl_window                   = config.bl_window
    e_grid_rt_trap              = config.e_grid_rt_trap
    enc_pickoff_trap            = config.enc_pickoff_trap

    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs, leftendpoint(bl_window), rightendpoint(bl_window))

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)

    # deconvolute waveform
    deconv_flt = InvCRFilter(τ)
    wvfs = deconv_flt.(wvfs)

    # get signal estimator
    signal_estimator = SignalEstimator(PolynomialDNI(config.kwargs_pars.sig_interpolation_order, config.kwargs_pars.sig_interpolation_length))

    # get energy grid for efficient optimization
    enc_trap_grid = zeros(Float64, length(e_grid_rt_trap), length(wvfs))
    for (r, rt) in enumerate(e_grid_rt_trap)
        if rt < ft
            continue
        end
        uflt_rtft      = TrapezoidalChargeFilter(rt, ft)
        
        enc_rtft       = signal_estimator.(uflt_rtft.(wvfs), enc_pickoff_trap)

        enc_trap_grid[r, :]   = enc_rtft
    end

    return enc_trap_grid
end
export dsp_trap_rt_optimization


"""
    dsp_cusp_rt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T},; ft::Quantity{T}=4.0u"µs") where T<:Real

Get ENC noise grid values for given CUSP grid rise times.

# Returns
    - `enc_cusp_grid`: Array ENC noise values for the given CUSP rise time grid
"""
function dsp_cusp_rt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T},; ft::Quantity{T}=2.0u"µs") where T<:Real
    # get config parameters
    bl_window                   = config.bl_window
    e_grid_rt_cusp              = config.e_grid_rt_cusp
    enc_pickoff_cusp            = config.enc_pickoff_cusp
    flt_length_cusp             = config.flt_length_cusp
    cusp_scale                  = ustrip(NoUnits, flt_length_cusp/step(wvfs[1].time))

    # set tau for CUSP filter to very high number to switch of CR filter
    τ_cusp = 10000000.0u"µs"

    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs, leftendpoint(bl_window), rightendpoint(bl_window))

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)

    # deconvolute waveform
    deconv_flt = InvCRFilter(τ)
    wvfs = deconv_flt.(wvfs)

    # get signal estimator
    signal_estimator = SignalEstimator(PolynomialDNI(config.kwargs_pars.sig_interpolation_order, config.kwargs_pars.sig_interpolation_length))

    # get energy grid for efficient optimization
    enc_cusp_grid = zeros(Float64, length(collect(e_grid_rt_cusp)), length(wvfs))
    for (r, rt) in enumerate(e_grid_rt_cusp)
        if rt < ft
            continue
        end
        uflt_rtft      = CUSPChargeFilter(rt, ft, τ_cusp, flt_length_cusp, cusp_scale)
        
        enc_rtft       = signal_estimator.(uflt_rtft.(wvfs), enc_pickoff_cusp)

        enc_cusp_grid[r, :]   = enc_rtft
    end

    return enc_cusp_grid
end
export dsp_cusp_rt_optimization


"""
    dsp_zac_rt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T},; ft::Quantity{T}=4.0u"µs") where T<:Real

Get ENC noise grid values for given ZAC grid rise times.

# Returns
    - `enc_zac_grid`: Array ENC noise values for the given ZAC rise time grid
"""
function dsp_zac_rt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T},; ft::Quantity{T}=2.0u"µs") where T<:Real
    # get config parameters
    bl_window                   = config.bl_window
    e_grid_rt_zac               = config.e_grid_rt_zac
    enc_pickoff_zac             = config.enc_pickoff_zac
    flt_length_zac              = config.flt_length_zac
    zac_scale                   = ustrip(NoUnits, flt_length_zac/step(wvfs[1].time))

    # set tau for ZAC filter to very high number to switch of CR filter
    τ_zac = 10000000.0u"µs"

    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs, leftendpoint(bl_window), rightendpoint(bl_window))

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)

    # deconvolute waveform
    deconv_flt = InvCRFilter(τ)
    wvfs = deconv_flt.(wvfs)

    # get signal estimator
    signal_estimator = SignalEstimator(PolynomialDNI(config.kwargs_pars.sig_interpolation_order, config.kwargs_pars.sig_interpolation_length))

    # get energy grid for efficient optimization
    enc_zac_grid = zeros(Float64, length(e_grid_rt_zac), length(wvfs))
    for (r, rt) in enumerate(e_grid_rt_zac)
        if rt < ft
            continue
        end
        uflt_rtft      = ZACChargeFilter(rt, ft, τ_zac, flt_length_zac, zac_scale)
        
        enc_rtft       = signal_estimator.(uflt_rtft.(wvfs), enc_pickoff_zac)

        enc_zac_grid[r, :]   = enc_rtft
    end

    return enc_zac_grid
end
export dsp_zac_rt_optimization


"""
    dsp_trap_ft_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, rt::Quantity{T}) where T<:Real

Get energy grid values for given trap grid rise times while varying the flat-top time.

# Returns
    - `e_grid`: Array energy values for the given trap rise time grid at a given flat-top time grid.
"""
function dsp_trap_ft_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, rt::Quantity{T}) where T<:Real
    # get config parameters
    bl_window                   = config.bl_window
    e_grid_ft_trap              = config.e_grid_ft_trap

    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs, leftendpoint(bl_window), rightendpoint(bl_window))

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)

    # deconvolute waveform
    deconv_flt = InvCRFilter(τ)
    wvfs = deconv_flt.(wvfs)

    # get signal estimator
    signal_estimator = SignalEstimator(PolynomialDNI(config.kwargs_pars.sig_interpolation_order, config.kwargs_pars.sig_interpolation_length))

    # t50 determination
    t50 = get_threshold(wvfs, maximum.(wvfs.signal) .* 0.5; mintot=config.kwargs_pars.tx_mintot)

    # get energy grid for efficient optimization
    e_grid   = Array{Union{Missing, Float32}}(missing, length(e_grid_ft_trap), length(wvfs))
    for (f, ft) in enumerate(e_grid_ft_trap)
        if rt < ft
            continue
        end
        uflt_rtft      = TrapezoidalChargeFilter(rt, ft)
        
        e_rtft         = signal_estimator.(uflt_rtft.(wvfs), t50 .+ (rt + ft/2))

        e_grid[f, :]     = e_rtft
    end

    return e_grid
end
export dsp_trap_ft_optimization


"""
    dsp_cusp_ft_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, rt::Quantity{T}) where T<:Real

Get energy grid values for given CUSP grid rise times while varying the flat-top time.

# Returns
    - `e_grid`: Array energy values for the given CUSP rise time grid at a given flat-top time grid.
"""
function dsp_cusp_ft_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, rt::Quantity{T}) where T<:Real
    # get config parameters
    bl_window                   = config.bl_window
    e_grid_ft_cusp              = config.e_grid_ft_cusp
    flt_length_cusp             = config.flt_length_cusp
    cusp_scale                  = ustrip(NoUnits, flt_length_cusp/step(wvfs[1].time))

    # set tau for CUSP filter to very high number to switch of CR filter
    τ_cusp = 10000000.0u"µs"

    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs, leftendpoint(bl_window), rightendpoint(bl_window))

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)

    # deconvolute waveform
    deconv_flt = InvCRFilter(τ)
    wvfs = deconv_flt.(wvfs)

    # get signal estimator
    signal_estimator = SignalEstimator(PolynomialDNI(config.kwargs_pars.sig_interpolation_order, config.kwargs_pars.sig_interpolation_length))

    # t50 determination
    t50 = get_threshold(wvfs, maximum.(wvfs.signal) .* 0.5; mintot=config.kwargs_pars.tx_mintot)

    # get energy grid for efficient optimization
    e_grid   = Array{Union{Missing, Float32}}(missing, length(e_grid_ft_cusp), length(wvfs))
    for (f, ft) in enumerate(e_grid_ft_cusp)
        if rt < ft
            continue
        end
        uflt_rtft      = CUSPChargeFilter(rt, ft, τ_cusp, flt_length_cusp, cusp_scale)
        
        e_rtft         = signal_estimator.(uflt_rtft.(wvfs), t50 .+ flt_length_cusp/2)

        e_grid[f, :]     = e_rtft
    end

    return e_grid
end
export dsp_cusp_ft_optimization


"""
    dsp_zac_ft_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, rt::Quantity{T}) where T<:Real

Get energy grid values for given ZAC grid rise times while varying the flat-top time.

Returns:
    - `e_grid`: Array energy values for the given ZAC rise time grid at a given flat-top time grid.
"""
function dsp_zac_ft_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, rt::Quantity{T}) where T<:Real
    # get config parameters
    bl_window                   = config.bl_window
    e_grid_ft_zac               = config.e_grid_ft_zac
    flt_length_zac              = config.flt_length_zac
    zac_scale                   = ustrip(NoUnits, flt_length_zac/step(wvfs[1].time))

    # set tau for ZAC filter to very high number to switch of CR filter
    τ_zac = 10000000.0u"µs"

    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs, leftendpoint(bl_window), rightendpoint(bl_window))

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)

    # deconvolute waveform
    deconv_flt = InvCRFilter(τ)
    wvfs = deconv_flt.(wvfs)

    # get signal estimator
    signal_estimator = SignalEstimator(PolynomialDNI(config.kwargs_pars.sig_interpolation_order, config.kwargs_pars.sig_interpolation_length))

    # t50 determination
    t50 = get_threshold(wvfs, maximum.(wvfs.signal) .* 0.5; mintot=config.kwargs_pars.tx_mintot)

    # get energy grid for efficient optimization
    e_grid   = Array{Union{Missing, Float32}}(missing, length(e_grid_ft_zac), length(wvfs))
    for (f, ft) in enumerate(e_grid_ft_zac)
        if rt < ft
            continue
        end
        uflt_rtft      = ZACChargeFilter(rt, ft, τ_zac, flt_length_zac, zac_scale)
        
        e_rtft         = signal_estimator.(uflt_rtft.(wvfs), t50 .+ flt_length_zac/2)

        e_grid[f, :]     = e_rtft
    end

    return e_grid
end
export dsp_zac_ft_optimization


""" 
    dsp_sg_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, pars_filter::PropDict) where T<:Real

Optimize the Savitzky-Golay filter parameters for a given waveform set.

# Returns
    - `aoe`: Array of efficiency values for the given Savitzky-Golay filter parameters
    - `e`: Array of energy values for the given Savitzky-Golay filter parameters
    - `blmean`: Baseline mean value
    - `blslope`: Baseline slope value
"""
function dsp_sg_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, pars_filter::PropDict) where T<:Real
    # get config parameters
    bl_window      = config.bl_window
    a_grid_wl_sg   = config.a_grid_wl_sg
    sg_flt_degree  = config.sg_flt_degree
    current_window = config.current_window

    # get optimal filter parameters
    rt = pars_filter.trap.rt
    ft = pars_filter.trap.ft

    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs, leftendpoint(bl_window), rightendpoint(bl_window))

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)

    # deconvolute waveform
    deconv_flt = InvCRFilter(τ)
    wvfs = deconv_flt.(wvfs)
    
    # get signal estimator
    signal_estimator = SignalEstimator(PolynomialDNI(config.kwargs_pars.sig_interpolation_order, config.kwargs_pars.sig_interpolation_length))

    # t50 determination
    t50 = get_threshold(wvfs, maximum.(wvfs.signal) .* 0.5; mintot=config.kwargs_pars.tx_mintot)

    # get energy for filter parameters
    uflt_rtft = TrapezoidalChargeFilter(rt, ft)
    e_rtft    = signal_estimator.(uflt_rtft.(wvfs), t50 .+ (rt + ft/2))

    # extract current with filter length in grid with second order polynominal and first derivative
    aoe_grid   = ones(Float64, length(a_grid_wl_sg), length(wvfs))
    for (w, wl) in enumerate(a_grid_wl_sg)
        # extract current with optimal SG filter length with second order polynominal and first derivative
        wvfs_sgflt_deriv = SavitzkyGolayFilter(wl, sg_flt_degree, 1).(wvfs)
        current_max = get_wvf_maximum.(wvfs_sgflt_deriv, leftendpoint(current_window), rightendpoint(current_window))

        aoe_grid[w, :]     = ustrip.(current_max) ./ e_rtft
    end
    return (aoe = aoe_grid, e = e_rtft, blmean = bl_stats.mean, blslope = bl_stats.slope, t50 = t50)
end
export dsp_sg_optimization


""" 
    dsp_sg_optimization_compressed(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, pars_filter::PropDict) where T<:Real

Optimize the Savitzky-Golay filter parameters for a given waveform set.

# Returns
    - `aoe`: Array of efficiency values for the given Savitzky-Golay filter parameters
    - `e`: Array of energy values for the given Savitzky-Golay filter parameters
    - `blmean`: Baseline mean value
    - `blslope`: Baseline slope value
"""
function dsp_sg_optimization_compressed(wvfs_wdw::ArrayOfRDWaveforms, wvfs_pre::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, pars_filter::PropDict; presum_rate::T = T(8)) where T<:Real
    # get config parameters
    bl_window      = config.bl_window
    a_grid_wl_sg   = config.a_grid_wl_sg
    sg_flt_degree  = config.sg_flt_degree
    current_window = config.current_window

    # get optimal filter parameters
    rt = pars_filter.trap.rt
    ft = pars_filter.trap.ft

    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs_pre, leftendpoint(bl_window), rightendpoint(bl_window))

    # substract baseline from waveforms
    wvfs_pre = shift_waveform.(wvfs_pre, -bl_stats.mean)
    wvfs_wdw = shift_waveform.(wvfs_wdw, -bl_stats.mean ./ presum_rate)

    # deconvolute waveform
    deconv_flt = InvCRFilter(τ)
    wvfs_pre = deconv_flt.(wvfs_pre)
    wvfs_wdw = deconv_flt.(wvfs_wdw) 

    # get wvf maximum
    wvf_max_pre = maximum.(wvfs_pre.signal)
    
    # get signal estimator
    signal_estimator = SignalEstimator(PolynomialDNI(config.kwargs_pars.sig_interpolation_order, config.kwargs_pars.sig_interpolation_length))

    # t50 determination
    t50_pre = get_threshold(wvfs_pre, wvf_max_pre .* 0.5; mintot=config.kwargs_pars.tx_mintot)

    # get energy for filter parameters
    uflt_rtft = TrapezoidalChargeFilter(rt, ft)
    e_rtft    = signal_estimator.(uflt_rtft.(wvfs_pre), t50_pre .+ (rt + ft/2))

    # extract current with filter length in grid with second order polynominal and first derivative
    aoe_grid   = ones(Float64, length(a_grid_wl_sg), length(wvfs_wdw))
    for (w, wl) in enumerate(a_grid_wl_sg)
        # extract current with optimal SG filter length with second order polynominal and first derivative
        a_sg = get_wvf_maximum.(SavitzkyGolayFilter(wl, sg_flt_degree, 1).(wvfs_wdw), leftendpoint(current_window), rightendpoint(current_window))
        aoe_grid[w, :]     = ustrip.(a_sg) ./ e_rtft
    end
    return (aoe = aoe_grid, e = e_rtft, blmean = bl_stats.mean, blslope = bl_stats.slope, t50 = t50_pre)
end
export dsp_sg_optimization_compressed