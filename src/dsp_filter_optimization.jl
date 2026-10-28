# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    dsp_qc_flt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Unitful.Time{<:Real}, f_evaluate_qc::Union{Function, Missing})

Get QC DSP for filter parameter optimization for a given waveform set.
If f_evaluate_qc is missing, skips expensive ML classification and sets qc_label=-1.
"""
function dsp_qc_flt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Unitful.Time{<:Real}, f_evaluate_qc::Function)
    _get_dsp_qc_flt_optimization(wvfs, config, τ, Base.Fix2(get_qc_classifier, f_evaluate_qc))
end
function dsp_qc_flt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Unitful.Time{<:Real}, ::Missing)
    _get_dsp_qc_flt_optimization(wvfs, config, τ, nothing)
end
export dsp_qc_flt_optimization

"""
    dsp_qc_flt_optimization_compressed(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Unitful.Time{<:Real}, f_evaluate_qc::Union{Function, Missing})

Get QC DSP for filter parameter optimization for a given waveform set.
If f_evaluate_qc is missing, skips expensive ML classification and sets qc_label=-1.
"""
function dsp_qc_flt_optimization_compressed(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Unitful.Time{<:Real}, f_evaluate_qc::Function) 
    _get_dsp_qc_flt_optimization(wvfs, config, τ, Base.Fix2(get_qc_classifier_compressed, f_evaluate_qc))
end
function dsp_qc_flt_optimization_compressed(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Unitful.Time{<:Real}, ::Missing) 
    _get_dsp_qc_flt_optimization(wvfs, config, τ, nothing)
end
export dsp_qc_flt_optimization_compressed

function _get_dsp_qc_flt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Unitful.Time{<:Real}, get_qc::Union{Function, Nothing})
    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs, leftendpoint(config.bl_window), rightendpoint(config.bl_window))

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)

    # deconvolute waveform
    deconv_flt = InvCRFilter(τ)
    wvfs = deconv_flt.(wvfs)

    # get signal estimator
    signal_estimator = SignalEstimator(PolynomialDNI(config.kwargs_pars.sig_interpolation_order, config.kwargs_pars.sig_interpolation_length))

    # t50 determination
    t50 = get_threshold(wvfs, maximum.(wvfs.signal) .* 0.5; mintot=config.kwargs_pars.tx_mintot)

    # get QC classifier labels (skip expensive Haar filtering if no ML model)
    qc_labels = isnothing(get_qc) ? fill(-1, length(wvfs)) : Int.(get_qc(wvfs))

    # get default filter parameters
    rt = config.default_flt_param.trap.rt
    ft = config.default_flt_param.trap.ft

    # get energy for simple trap filtering
    uflt_rtft = TrapezoidalChargeFilter(rt, ft)
    e_rtft    = signal_estimator.(uflt_rtft.(wvfs), t50 .+ (rt + ft/2))

    return TypedTables.Table(
                energy = e_rtft, 
                blmean = bl_stats.mean, blslope = bl_stats.slope, 
                t50 = t50,
                qc_label = qc_labels
            )
end

"""
    dsp_qdrift_flt_optimization(wvfs::ArrayOfRDWaveforms, blmean::Vector{<:Real}, config::DSPConfig, τ::Unitful.Time{<:Real})

Get QDrift filter parameter optimization for a given waveform set.
"""
function dsp_qdrift_flt_optimization(wvfs::ArrayOfRDWaveforms, blmean::Vector{<:Real}, config::DSPConfig, τ::Unitful.Time{<:Real})
    # get config parameters
    t0_threshold      = config.t0_threshold
    qdrift_int_length = config.qdrift_int_length

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -blmean)

    # deconvolute waveform 
    # --> wvfs = wvfs_pz
    deconv_flt = InvCRFilter(τ)
    wvfs = deconv_flt.(wvfs)

    # t0 determination
    t0 = get_t0(wvfs, t0_threshold; flt_pars=config.kwargs_pars.t0_flt_pars, mintot=config.kwargs_pars.t0_mintot)
    
    # get qdrift
    get_qdrift(wvfs, t0, qdrift_int_length; pol_power=config.kwargs_pars.int_interpolation_order, sign_est_length=config.kwargs_pars.int_interpolation_length)
end
export dsp_qdrift_flt_optimization


"""
    dsp_flt_optimization_prepare(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Unitful.Time{<:Real})

Baseline-subtract and pole-zero correct `wvfs` once for repeated filter evaluations.

# Returns
A `NamedTuple` with
- `wvfs_pz`: baseline-subtracted waveforms deconvolved with `InvCRFilter(τ)`
- `t50`: time of 50% of the maximum of each `wvfs_pz`
- `bl_stats`: baseline statistics over `config.bl_window`
- `signal_estimator`: `SignalEstimator` used to read the filter output at the pickoff time
- `sampling_period`: time step of the waveforms
"""
function dsp_flt_optimization_prepare(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Unitful.Time{<:Real})
    bl_stats = signalstats.(wvfs, leftendpoint(config.bl_window), rightendpoint(config.bl_window))
    wvfs_pz = InvCRFilter(τ).(shift_waveform.(wvfs, -bl_stats.mean))
    t50 = get_threshold(wvfs_pz, maximum.(wvfs_pz.signal) .* 0.5; mintot=config.kwargs_pars.tx_mintot)
    signal_estimator = SignalEstimator(PolynomialDNI(config.kwargs_pars.sig_interpolation_order, config.kwargs_pars.sig_interpolation_length))
    return (; wvfs_pz, t50, bl_stats, signal_estimator, sampling_period=step(first(wvfs).time))
end
export dsp_flt_optimization_prepare


"""
    dsp_energy_rtft(prep::NamedTuple, config::DSPConfig, rt, ft, flt)

Energy of the prepared waveforms `prep` (see [`dsp_flt_optimization_prepare`](@ref))
at rise time `rt` and flat-top time `ft`. `flt` is `Val(:trap)`, `Val(:cusp)`, `Val(:zac)`
or the corresponding `Symbol`.

The trapezoidal energy is read at `t50 + rt + ft/2`, the CUSP and ZAC energies at
`t50 + flt_length/2`. CUSP and ZAC use `τ = 10000000 µs`, which disables their CR filter,
and the filter length of `config`.

Each call allocates transient memory of about five times the size of `prep.wvfs_pz` for
CUSP and ZAC and about twice for the trapezoidal filter; callers running it concurrently
must bound the number of simultaneous calls accordingly. `rt` and `ft` may be given in any
time unit and number type; they are converted to the type of the sampling period (trapezoidal)
or of the filter length (CUSP, ZAC), in floating point, before the filter is built. An energy is `NaN` for an event
whose pickoff time lies outside the filtered time axis. The function only reads `prep` and
is safe to call concurrently on the same `prep`.
"""
function dsp_energy_rtft(prep::NamedTuple, config::DSPConfig, rt::Unitful.Time{<:Real}, ft::Unitful.Time{<:Real}, ::Val{:trap})
    T = typeof(float(prep.sampling_period))
    uflt_rtft = TrapezoidalChargeFilter(convert(T, rt), convert(T, ft))
    return prep.signal_estimator.(uflt_rtft.(prep.wvfs_pz), prep.t50 .+ (rt + ft/2))
end

function dsp_energy_rtft(prep::NamedTuple, config::DSPConfig, rt::Unitful.Time{<:Real}, ft::Unitful.Time{<:Real}, ::Val{:cusp})
    flt_length = config.flt_length_cusp
    scale = ustrip(NoUnits, flt_length / prep.sampling_period)
    T = typeof(float(flt_length))
    uflt_rtft = CUSPChargeFilter(convert(T, rt), convert(T, ft), convert(T, 10000000.0u"µs"), flt_length, scale)
    return prep.signal_estimator.(uflt_rtft.(prep.wvfs_pz), prep.t50 .+ flt_length/2)
end

function dsp_energy_rtft(prep::NamedTuple, config::DSPConfig, rt::Unitful.Time{<:Real}, ft::Unitful.Time{<:Real}, ::Val{:zac})
    flt_length = config.flt_length_zac
    scale = ustrip(NoUnits, flt_length / prep.sampling_period)
    T = typeof(float(flt_length))
    uflt_rtft = ZACChargeFilter(convert(T, rt), convert(T, ft), convert(T, 10000000.0u"µs"), flt_length, scale)
    return prep.signal_estimator.(uflt_rtft.(prep.wvfs_pz), prep.t50 .+ flt_length/2)
end

dsp_energy_rtft(prep::NamedTuple, config::DSPConfig, rt::Unitful.Time{<:Real}, ft::Unitful.Time{<:Real}, flt::Symbol) =
    dsp_energy_rtft(prep, config, rt, ft, Val(flt))
export dsp_energy_rtft

"""
    realized_flattop(ft, flt_length, sampling_period, flt)

Flat-top time that the filter `flt` (`Val(:trap)`, `Val(:cusp)`, `Val(:zac)` or the
corresponding `Symbol`) applies when it is built with the flat-top time `ft` for waveforms
sampled with `sampling_period`. `flt_length` is the filter length of CUSP and ZAC and is not
read for the trapezoidal filter.

`RadiationDetectorDSP` rounds the times of the filters to whole samples:
- trapezoidal filter: the gap is `round(ft / sampling_period)` samples;
- CUSP and ZAC: the flat top is `FT = round(ft / sampling_period)` samples and the filter length
  is `N = round(flt_length / sampling_period)` samples. The two sides of the filter have
  `(N - FT) / 2` samples each, so the flat top is lengthened by one sample if `N - FT` is odd.

The rule mirrors the installed `RadiationDetectorDSP` (`cusp_charge_filter_coeffs`,
`zac_charge_filter_coeffs`, `fltinstance`). Two flat-top times with the same realized value
produce identical filters, and the realized value is its own realized value.

Returns the realized flat-top time, a multiple of `sampling_period`.
"""
function realized_flattop(ft::Unitful.Time{<:Real}, flt_length::Unitful.Time{<:Real}, sampling_period::Unitful.Time{<:Real}, ::Val{:trap})
    return round(Int, ustrip(NoUnits, ft / sampling_period)) * sampling_period
end

function realized_flattop(ft::Unitful.Time{<:Real}, flt_length::Unitful.Time{<:Real}, sampling_period::Unitful.Time{<:Real}, ::Union{Val{:cusp}, Val{:zac}})
    n_ft = round(Int, ustrip(NoUnits, ft / sampling_period))
    n_flt = round(Int, ustrip(NoUnits, flt_length / sampling_period))
    return (iseven(n_flt - n_ft) ? n_ft : n_ft + 1) * sampling_period
end

realized_flattop(ft::Unitful.Time{<:Real}, flt_length::Unitful.Time{<:Real}, sampling_period::Unitful.Time{<:Real}, flt::Symbol) =
    realized_flattop(ft, flt_length, sampling_period, Val(flt))
export realized_flattop

"""
    dsp_pulse_shape_flt_optimization(prep::NamedTuple, config::DSPConfig, sg_wl::Unitful.Time{<:Real})

Pulse-shape quantities of the baseline-subtracted, deconvolved waveforms `prep.wvfs_pz`; only
this field of `prep` is read. Pass waveforms whose sampling resolves the SG window `sg_wl`, as
the windowed waveforms of the dsp tier do (for example `(wvfs_pz = wvfs_wdw_pz,)`); the
presummed waveforms prepared by [`dsp_flt_optimization_prepare`](@ref) can be sampled more
coarsely than `sg_wl`, in which case `a` is zero.

# Returns
A `NamedTuple` with
- `t10`, `t90`: times of 10% and 90% of the maximum of the waveforms in `prep.wvfs_pz`
- `a`: maximum of the Savitzky-Golay current (window length `sg_wl`, polynomial degree
  `config.sg_flt_degree`) within `config.current_window`
"""
function dsp_pulse_shape_flt_optimization(prep::NamedTuple, config::DSPConfig, sg_wl::Unitful.Time{<:Real})
    wvf_max = maximum.(prep.wvfs_pz.signal)
    t10 = get_threshold(prep.wvfs_pz, wvf_max .* 0.1; mintot=config.kwargs_pars.tx_mintot)
    t90 = get_threshold(prep.wvfs_pz, wvf_max .* 0.9; mintot=config.kwargs_pars.tx_mintot)
    wvfs_sgflt_deriv = SavitzkyGolayFilter(sg_wl, config.sg_flt_degree, 1).(prep.wvfs_pz)
    a = get_wvf_maximum.(wvfs_sgflt_deriv, leftendpoint(config.current_window), rightendpoint(config.current_window)).max
    return (; t10, t90, a)
end
export dsp_pulse_shape_flt_optimization


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

    # set τ for CUSP filter to very high number to switch of CR filter
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

    # set τ for ZAC filter to very high number to switch of CR filter
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
    prep = dsp_flt_optimization_prepare(wvfs, config, τ)
    e_grid = Array{Union{Missing, Float32}}(missing, length(config.e_grid_ft_trap), length(wvfs))
    for (f, ft) in pairs(config.e_grid_ft_trap)
        e_grid[f, :] = dsp_energy_rtft(prep, config, rt, ft, Val(:trap))
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
    prep = dsp_flt_optimization_prepare(wvfs, config, τ)
    e_grid = Array{Union{Missing, Float32}}(missing, length(config.e_grid_ft_cusp), length(wvfs))
    for (f, ft) in pairs(config.e_grid_ft_cusp)
        e_grid[f, :] = dsp_energy_rtft(prep, config, rt, ft, Val(:cusp))
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
    prep = dsp_flt_optimization_prepare(wvfs, config, τ)
    e_grid = Array{Union{Missing, Float32}}(missing, length(config.e_grid_ft_zac), length(wvfs))
    for (f, ft) in pairs(config.e_grid_ft_zac)
        e_grid[f, :] = dsp_energy_rtft(prep, config, rt, ft, Val(:zac))
    end
    return e_grid
end
export dsp_zac_ft_optimization


""" 
    dsp_sg_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, pars_filter::PropDict; f_evaluate_qc::Union{Function, Missing}=missing) where T<:Real

Optimize the Savitzky-Golay filter parameters for a given waveform set.

# Arguments
    - `f_evaluate_qc`: Optional ML QC evaluation function. If missing, default qc_label=-1 is used.

# Returns
    - `aoe`: Array of efficiency values for the given Savitzky-Golay filter parameters
    - `energy`: Array of energy values for the given Savitzky-Golay filter parameters
    - `blmean`: Baseline mean value
    - `blslope`: Baseline slope value
    - `qc_label`: QC labels (-1 if no ML model available)
"""
function dsp_sg_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, pars_filter::PropDict; f_evaluate_qc::Union{Function, Missing}=missing) where T<:Real
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

    # get QC classifier labels (skip expensive ML classification if no model)
    qc_labels = !ismissing(f_evaluate_qc) ? Int.(get_qc_classifier(wvfs, f_evaluate_qc)) : fill(-1, length(wvfs))

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
        current_max = get_wvf_maximum.(wvfs_sgflt_deriv, leftendpoint(current_window), rightendpoint(current_window)).max

        aoe_grid[w, :]     = ustrip.(current_max) ./ e_rtft
    end
    return TypedTables.Table(aoe = VectorOfSimilarVectors(aoe_grid), energy = e_rtft, 
        blmean = bl_stats.mean, blslope = bl_stats.slope, 
        t50 = t50,
        qc_label = qc_labels
    )
end
export dsp_sg_optimization


""" 
    dsp_sg_optimization_compressed(wvfs_wdw::ArrayOfRDWaveforms, wvfs_pre::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, pars_filter::PropDict; presum_rate::Real=T(8), f_evaluate_qc::Union{Function, Missing}=missing) where T<:Real

Optimize the Savitzky-Golay filter parameters for a given waveform set.

# Arguments
    - `f_evaluate_qc`: Optional ML QC evaluation function. If missing, default qc_label=-1 is used.

# Returns
    - `aoe`: Array of efficiency values for the given Savitzky-Golay filter parameters
    - `energy`: Array of energy values for the given Savitzky-Golay filter parameters
    - `blmean`: Baseline mean value
    - `blslope`: Baseline slope value
    - `qc_label`: QC labels (-1 if no ML model available)
"""
function dsp_sg_optimization_compressed(wvfs_wdw::ArrayOfRDWaveforms, wvfs_pre::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, pars_filter::PropDict; presum_rate::Real=T(8), f_evaluate_qc::Union{Function, Missing}=missing) where T<:Real
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

    # get QC classifier labels (skip expensive ML classification if no model)
    qc_labels = !ismissing(f_evaluate_qc) ? Int.(get_qc_classifier_compressed(wvfs_pre, f_evaluate_qc)) : fill(-1, length(wvfs_pre))

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
        a_sg = get_wvf_maximum.(SavitzkyGolayFilter(wl, sg_flt_degree, 1).(wvfs_wdw), leftendpoint(current_window), rightendpoint(current_window)).max
        aoe_grid[w, :]     = ustrip.(a_sg) ./ e_rtft
    end
    return TypedTables.Table(aoe = VectorOfSimilarVectors(aoe_grid), energy = e_rtft, 
                blmean = bl_stats.mean, blslope = bl_stats.slope, 
                t50 = t50_pre,
                qc_label = qc_labels
                )
end
export dsp_sg_optimization_compressed

"""
    dsp_mwa_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, pars_filter::PropDict; f_evaluate_qc::Union{Function, Missing}=missing) where T<:Real

Optimize the window length of the MWA current filter (`MovingWindowMultiFilter` on the
simple-derivative current) for the A/E parameter, same scheme as `dsp_sg_optimization`.

# Arguments
    - `f_evaluate_qc`: Optional ML QC evaluation function. If missing, default qc_label=-1 is used.

# Returns
    - `aoe`: A/E values for each MWA window length in `config.a_grid_wl_mwa`
    - `energy`: Array of energy values
    - `blmean`: Baseline mean value
    - `blslope`: Baseline slope value
    - `t50`: 50% rise time of the waveform
    - `qc_label`: QC labels (-1 if no ML model available)
"""
function dsp_mwa_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, pars_filter::PropDict; f_evaluate_qc::Union{Function, Missing}=missing) where T<:Real
    # get config parameters
    bl_window      = config.bl_window
    a_grid_wl_mwa  = config.a_grid_wl_mwa
    current_window = config.current_window

    # get optimal filter parameters
    rt = pars_filter.trap.rt
    ft = pars_filter.trap.ft

    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs, leftendpoint(bl_window), rightendpoint(bl_window))

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)

    # get QC classifier labels (skip expensive ML classification if no model)
    qc_labels = !ismissing(f_evaluate_qc) ? Int.(get_qc_classifier(wvfs, f_evaluate_qc)) : fill(-1, length(wvfs))

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

    # current from the simple derivative, MWA-shaped with the window length in the grid
    # (same chain as the a_mwa column in `dsp_icpc`)
    wvfs_deriv = DerivativeFilter(1).(wvfs)
    aoe_grid   = ones(Float64, length(a_grid_wl_mwa), length(wvfs))
    for (w, wl) in enumerate(a_grid_wl_mwa)
        a_mwa = get_wvf_maximum.(MovingWindowMultiFilter(wl).(wvfs_deriv), leftendpoint(current_window), rightendpoint(current_window)).max
        aoe_grid[w, :] = ustrip.(a_mwa) ./ e_rtft
    end
    return TypedTables.Table(aoe = VectorOfSimilarVectors(aoe_grid), energy = e_rtft,
        blmean = bl_stats.mean, blslope = bl_stats.slope,
        t50 = t50,
        qc_label = qc_labels
    )
end
export dsp_mwa_optimization

"""
    dsp_mwa_optimization_compressed(wvfs_wdw::ArrayOfRDWaveforms, wvfs_pre::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, pars_filter::PropDict; presum_rate::Real=T(8), f_evaluate_qc::Union{Function, Missing}=missing)

Optimize the window length of the MWA current filter (`MovingWindowMultiFilter` on the
simple-derivative current) for the A/E parameter, same scheme as `dsp_sg_optimization_compressed`.

# Returns
    - `aoe`: A/E values for each MWA window length in `config.a_grid_wl_mwa`
    - `energy`: Array of energy values
    - `blmean`: Baseline mean value
    - `blslope`: Baseline slope value
    - `t50`: 50% rise time of the presummed waveform
    - `qc_label`: QC labels (-1 if no ML model available)
"""
function dsp_mwa_optimization_compressed(wvfs_wdw::ArrayOfRDWaveforms, wvfs_pre::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, pars_filter::PropDict; presum_rate::Real=T(8), f_evaluate_qc::Union{Function, Missing}=missing) where T<:Real
    # get config parameters
    bl_window      = config.bl_window
    a_grid_wl_mwa  = config.a_grid_wl_mwa
    current_window = config.current_window

    # get optimal filter parameters
    rt = pars_filter.trap.rt
    ft = pars_filter.trap.ft

    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs_pre, leftendpoint(bl_window), rightendpoint(bl_window))

    # substract baseline from waveforms
    wvfs_pre = shift_waveform.(wvfs_pre, -bl_stats.mean)
    wvfs_wdw = shift_waveform.(wvfs_wdw, -bl_stats.mean ./ presum_rate)

    # get QC classifier labels (skip expensive ML classification if no model)
    qc_labels = !ismissing(f_evaluate_qc) ? Int.(get_qc_classifier_compressed(wvfs_pre, f_evaluate_qc)) : fill(-1, length(wvfs_pre))

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

    # current from the simple derivative, MWA-shaped with the window length in the grid
    # (same chain as the a_mwa column in `dsp_icpc`)
    wvfs_deriv = DerivativeFilter(1).(wvfs_wdw)
    aoe_grid   = ones(Float64, length(a_grid_wl_mwa), length(wvfs_wdw))
    for (w, wl) in enumerate(a_grid_wl_mwa)
        a_mwa = get_wvf_maximum.(MovingWindowMultiFilter(wl).(wvfs_deriv), leftendpoint(current_window), rightendpoint(current_window)).max
        aoe_grid[w, :] = ustrip.(a_mwa) ./ e_rtft
    end
    return TypedTables.Table(aoe = VectorOfSimilarVectors(aoe_grid), energy = e_rtft,
                blmean = bl_stats.mean, blslope = bl_stats.slope,
                t50 = t50_pre,
                qc_label = qc_labels
                )
end
export dsp_mwa_optimization_compressed
