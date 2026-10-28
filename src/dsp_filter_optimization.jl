# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    dsp_qc_flt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Unitful.Time, f_evaluate_qc::Function)
    dsp_qc_flt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Unitful.Time, f_evaluate_qc::Missing)

Energy, baseline and quality-cut label of each waveform with the default trapezoidal filter, used to select waveforms for the filter-parameter scans.

The waveforms are decoded, the baseline mean in `config.bl_window` is subtracted, and the waveforms are pole-zero corrected with the decay time `τ`. The energy is the value of a `TrapezoidalChargeFilter(rt, ft)` with `rt` and `ft` of `config.default_flt_param.trap`, determined by a signal estimator at `t50 + rt + ft/2`. `t50` is the time at which the waveform reaches 50% of its maximum.

The label is `get_qc_classifier(wvfs, f_evaluate_qc)` of the baseline-subtracted and pole-zero-corrected waveforms. With `missing` as `f_evaluate_qc` no classification is done and `qc_label` is -1.

# Returns
- table with the columns `energy`, `blmean`, `blslope`, `t50` (in µs) and `qc_label`, one row per waveform
"""
function dsp_qc_flt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Unitful.Time{<:Real}, f_evaluate_qc::Function)
    _get_dsp_qc_flt_optimization(wvfs, config, τ, Base.Fix2(get_qc_classifier, f_evaluate_qc))
end
function dsp_qc_flt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Unitful.Time{<:Real}, ::Missing)
    _get_dsp_qc_flt_optimization(wvfs, config, τ, nothing)
end
export dsp_qc_flt_optimization

"""
    dsp_qc_flt_optimization_compressed(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Unitful.Time, f_evaluate_qc::Function)
    dsp_qc_flt_optimization_compressed(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Unitful.Time, f_evaluate_qc::Missing)

Variant of [`dsp_qc_flt_optimization`](@ref) for presummed waveforms. The label is determined with [`get_qc_classifier_compressed`](@ref); everything else is the same.
"""
function dsp_qc_flt_optimization_compressed(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Unitful.Time{<:Real}, f_evaluate_qc::Function) 
    _get_dsp_qc_flt_optimization(wvfs, config, τ, Base.Fix2(get_qc_classifier_compressed, f_evaluate_qc))
end
function dsp_qc_flt_optimization_compressed(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Unitful.Time{<:Real}, ::Missing) 
    _get_dsp_qc_flt_optimization(wvfs, config, τ, nothing)
end
export dsp_qc_flt_optimization_compressed

function _get_dsp_qc_flt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Unitful.Time{<:Real}, get_qc::Union{Function, Nothing})
    # get waveform data
    wvfs = decode_data(wvfs)

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
    dsp_qdrift_flt_optimization(wvfs::ArrayOfRDWaveforms, blmean::Vector{<:Real}, config::DSPConfig, τ::Unitful.Time)

Q-drift parameter of each waveform, for scanning the integration times `config.qdrift_int_length`.

The waveforms are decoded, `blmean` is subtracted, and the waveforms are pole-zero corrected with `τ`. Then t0 is determined with [`get_t0`](@ref) and the parameter with [`get_qdrift`](@ref).

# Arguments
- `blmean`: baseline mean of each waveform

# Returns
- vector of Q-drift values, one per waveform
"""
function dsp_qdrift_flt_optimization(wvfs::ArrayOfRDWaveforms, blmean::Vector{<:Real}, config::DSPConfig, τ::Unitful.Time{<:Real})
    # get config parameters
    t0_threshold      = config.t0_threshold
    qdrift_int_length = config.qdrift_int_length

    # get waveform data
    wvfs = decode_data(wvfs)

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
    dsp_trap_rt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}; ft::Quantity{T}=2.0u"µs") where T<:Real

Values of the trapezoidal filter at a fixed time for each rise time in `config.e_grid_rt_trap`, for estimating the equivalent noise charge (ENC).

The waveforms are decoded, the baseline mean in `config.bl_window` is subtracted, and the waveforms are pole-zero corrected with the decay time `τ`. For each rise time `rt` the waveforms are filtered with a `TrapezoidalChargeFilter(rt, ft)` and the filtered waveform is evaluated at `config.enc_pickoff_trap` by a signal estimator.

# Keywords
- `ft`: flat-top time

# Returns
- `Matrix{Float64}` of size `(length(config.e_grid_rt_trap), length(wvfs))`: row `r` holds the values for the `r`-th rise time
"""
function dsp_trap_rt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T},; ft::Quantity{T}=2.0u"µs") where T<:Real
    # get config parameters
    bl_window                   = config.bl_window
    e_grid_rt_trap              = config.e_grid_rt_trap
    enc_pickoff_trap            = config.enc_pickoff_trap

    # get waveform data
    wvfs = decode_data(wvfs)

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
    dsp_cusp_rt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}; ft::Quantity{T}=2.0u"µs") where T<:Real

Values of the CUSP filter at a fixed time for each rise time in `config.e_grid_rt_cusp`, for estimating the equivalent noise charge (ENC).

The waveforms are decoded, the baseline mean in `config.bl_window` is subtracted, and the waveforms are pole-zero corrected with the decay time `τ`. For each rise time `rt` the waveforms are filtered with a `CUSPChargeFilter` with the flat-top time `ft`, the length `config.flt_length_cusp` and the CR filter switched off and the filtered waveform is evaluated at `config.enc_pickoff_cusp` by a signal estimator.

# Keywords
- `ft`: flat-top time

# Returns
- `Matrix{Float64}` of size `(length(config.e_grid_rt_cusp), length(wvfs))`: row `r` holds the values for the `r`-th rise time
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

    # get waveform data
    wvfs = decode_data(wvfs)

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
    dsp_zac_rt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}; ft::Quantity{T}=2.0u"µs") where T<:Real

Values of the ZAC filter at a fixed time for each rise time in `config.e_grid_rt_zac`, for estimating the equivalent noise charge (ENC).

The waveforms are decoded, the baseline mean in `config.bl_window` is subtracted, and the waveforms are pole-zero corrected with the decay time `τ`. For each rise time `rt` the waveforms are filtered with a `ZACChargeFilter` with the flat-top time `ft`, the length `config.flt_length_zac` and the CR filter switched off and the filtered waveform is evaluated at `config.enc_pickoff_zac` by a signal estimator.

# Keywords
- `ft`: flat-top time

# Returns
- `Matrix{Float64}` of size `(length(config.e_grid_rt_zac), length(wvfs))`: row `r` holds the values for the `r`-th rise time
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

    # get waveform data
    wvfs = decode_data(wvfs)

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

Energy of the trapezoidal filter for the rise time `rt` and each flat-top time in `config.e_grid_ft_trap`.

The waveforms are decoded, the baseline mean in `config.bl_window` is subtracted, and the waveforms are pole-zero corrected with the decay time `τ`. `t50` is the time at which a waveform reaches 50% of its maximum. For each flat-top time `ft` the waveforms are filtered with a `TrapezoidalChargeFilter(rt, ft)`, and the energy is the value of the filtered waveform at `t50 + rt + ft/2` determined by a signal estimator.

# Returns
- `Matrix{Union{Missing, Float32}}` of size `(length(config.e_grid_ft_trap), length(wvfs))`: row `f` holds the energies for the `f`-th flat-top time
"""
function dsp_trap_ft_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, rt::Quantity{T}) where T<:Real
    # get config parameters
    bl_window                   = config.bl_window
    e_grid_ft_trap              = config.e_grid_ft_trap

    # get waveform data
    wvfs = decode_data(wvfs)

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
        
        uflt_rtft      = TrapezoidalChargeFilter(rt, ft)
        
        e_rtft         = signal_estimator.(uflt_rtft.(wvfs), t50 .+ (rt + ft/2))

        e_grid[f, :]     = e_rtft
    end

    return e_grid
end
export dsp_trap_ft_optimization


"""
    dsp_cusp_ft_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, rt::Quantity{T}) where T<:Real

Energy of the CUSP filter for the rise time `rt` and each flat-top time in `config.e_grid_ft_cusp`.

The waveforms are decoded, the baseline mean in `config.bl_window` is subtracted, and the waveforms are pole-zero corrected with the decay time `τ`. `t50` is the time at which a waveform reaches 50% of its maximum. For each flat-top time `ft` the waveforms are filtered with a `CUSPChargeFilter(rt, ft, ...)` of length `config.flt_length_cusp` with the CR filter switched off, and the energy is the value of the filtered waveform at `t50 + config.flt_length_cusp/2` determined by a signal estimator.

# Returns
- `Matrix{Union{Missing, Float32}}` of size `(length(config.e_grid_ft_cusp), length(wvfs))`: row `f` holds the energies for the `f`-th flat-top time
"""
function dsp_cusp_ft_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, rt::Quantity{T}) where T<:Real
    # get config parameters
    bl_window                   = config.bl_window
    e_grid_ft_cusp              = config.e_grid_ft_cusp
    flt_length_cusp             = config.flt_length_cusp
    cusp_scale                  = ustrip(NoUnits, flt_length_cusp/step(wvfs[1].time))

    # set τ for CUSP filter to very high number to switch of CR filter
    τ_cusp = 10000000.0u"µs"

    # get waveform data
    wvfs = decode_data(wvfs)

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
        
        uflt_rtft      = CUSPChargeFilter(rt, ft, τ_cusp, flt_length_cusp, cusp_scale)
        
        e_rtft         = signal_estimator.(uflt_rtft.(wvfs), t50 .+ flt_length_cusp/2)

        e_grid[f, :]     = e_rtft
    end

    return e_grid
end
export dsp_cusp_ft_optimization


"""
    dsp_zac_ft_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, rt::Quantity{T}) where T<:Real

Energy of the ZAC filter for the rise time `rt` and each flat-top time in `config.e_grid_ft_zac`.

The waveforms are decoded, the baseline mean in `config.bl_window` is subtracted, and the waveforms are pole-zero corrected with the decay time `τ`. `t50` is the time at which a waveform reaches 50% of its maximum. For each flat-top time `ft` the waveforms are filtered with a `ZACChargeFilter(rt, ft, ...)` of length `config.flt_length_zac` with the CR filter switched off, and the energy is the value of the filtered waveform at `t50 + config.flt_length_zac/2` determined by a signal estimator.

# Returns
- `Matrix{Union{Missing, Float32}}` of size `(length(config.e_grid_ft_zac), length(wvfs))`: row `f` holds the energies for the `f`-th flat-top time
"""
function dsp_zac_ft_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, rt::Quantity{T}) where T<:Real
    # get config parameters
    bl_window                   = config.bl_window
    e_grid_ft_zac               = config.e_grid_ft_zac
    flt_length_zac              = config.flt_length_zac
    zac_scale                   = ustrip(NoUnits, flt_length_zac/step(wvfs[1].time))

    # set τ for ZAC filter to very high number to switch of CR filter
    τ_zac = 10000000.0u"µs"

    # get waveform data
    wvfs = decode_data(wvfs)

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
        
        uflt_rtft      = ZACChargeFilter(rt, ft, τ_zac, flt_length_zac, zac_scale)
        
        e_rtft         = signal_estimator.(uflt_rtft.(wvfs), t50 .+ flt_length_zac/2)

        e_grid[f, :]     = e_rtft
    end

    return e_grid
end
export dsp_zac_ft_optimization


"""
    dsp_sg_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, pars_filter::PropDict; f_evaluate_qc::Union{Function, Missing}=missing) where T<:Real

Energy and A/E for each window length in `config.a_grid_wl_sg` of the Savitzky-Golay current filter.

The waveforms are decoded and the baseline mean in `config.bl_window` is subtracted. The quality-cut label is determined from these waveforms. The waveforms are then pole-zero corrected with `τ`. The energy is the value of a `TrapezoidalChargeFilter(rt, ft)` at `t50 + rt + ft/2`, determined by a signal estimator, where `t50` is the time at which the waveform reaches 50% of its maximum. For each window length `wl` the current is the maximum in `config.current_window` of the waveform filtered with `SavitzkyGolayFilter(wl, config.sg_flt_degree, 1)`. A/E is the current without its unit divided by the energy.

# Arguments
- `pars_filter`: needs `pars_filter.trap.rt` and `pars_filter.trap.ft`, the rise and flat-top time of the energy filter

# Keywords
- `f_evaluate_qc`: function for [`get_qc_classifier`](@ref). With `missing` no classification is done and `qc_label` is -1.

# Returns
- table with one row per waveform and the columns `aoe` (A/E for all window lengths of the grid, as a vector per waveform), `energy`, `blmean`, `blslope`, `t50` (in µs) and `qc_label`
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

    # get waveform data
    wvfs = decode_data(wvfs)

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

Variant of [`dsp_sg_optimization`](@ref) for compressed data.

The baseline, quality-cut label, `t50` and energy are determined from the presummed waveforms `wvfs_pre`. The current is determined from the windowed waveforms `wvfs_wdw`, from which `blmean / presum_rate` is subtracted. The output table has the same columns; `t50` is the time of the presummed waveform.

# Keywords
- `presum_rate`: number of samples summed in the presummed waveforms
- `f_evaluate_qc`: function for [`get_qc_classifier_compressed`](@ref)
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

    # get waveform data
    wvfs_wdw = decode_data(wvfs_wdw)
    wvfs_pre = decode_data(wvfs_pre)

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

Energy and A/E for each window length in `config.a_grid_wl_mwa` of the moving-window-average current filter.

Same as [`dsp_sg_optimization`](@ref), except for the current: the pole-zero-corrected waveforms are filtered with `DerivativeFilter(1)` and then with `MovingWindowMultiFilter(wl; passes = config.a_mwa_passes)` for each window length `wl`, and the maximum in `config.current_window` is the current. This is the chain of the column `a_mwa` of [`dsp_icpc`](@ref).
"""
function dsp_mwa_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, pars_filter::PropDict; f_evaluate_qc::Union{Function, Missing}=missing) where T<:Real
    # get config parameters
    bl_window      = config.bl_window
    a_grid_wl_mwa  = config.a_grid_wl_mwa
    current_window = config.current_window

    # get optimal filter parameters
    rt = pars_filter.trap.rt
    ft = pars_filter.trap.ft

    # get waveform data
    wvfs = decode_data(wvfs)

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
        a_mwa = get_wvf_maximum.(MovingWindowMultiFilter(wl; passes = config.a_mwa_passes).(wvfs_deriv), leftendpoint(current_window), rightendpoint(current_window)).max
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
    dsp_mwa_optimization_compressed(wvfs_wdw::ArrayOfRDWaveforms, wvfs_pre::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, pars_filter::PropDict; presum_rate::Real=T(8), f_evaluate_qc::Union{Function, Missing}=missing) where T<:Real

Variant of [`dsp_mwa_optimization`](@ref) for compressed data.

The arguments and the split of the work between `wvfs_pre` and `wvfs_wdw` are the same as for [`dsp_sg_optimization_compressed`](@ref).
"""
function dsp_mwa_optimization_compressed(wvfs_wdw::ArrayOfRDWaveforms, wvfs_pre::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, pars_filter::PropDict; presum_rate::Real=T(8), f_evaluate_qc::Union{Function, Missing}=missing) where T<:Real
    # get config parameters
    bl_window      = config.bl_window
    a_grid_wl_mwa  = config.a_grid_wl_mwa
    current_window = config.current_window

    # get optimal filter parameters
    rt = pars_filter.trap.rt
    ft = pars_filter.trap.ft

    # get waveform data
    wvfs_wdw = decode_data(wvfs_wdw)
    wvfs_pre = decode_data(wvfs_pre)

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
        a_mwa = get_wvf_maximum.(MovingWindowMultiFilter(wl; passes = config.a_mwa_passes).(wvfs_deriv), leftendpoint(current_window), rightendpoint(current_window)).max
        aoe_grid[w, :] = ustrip.(a_mwa) ./ e_rtft
    end
    return TypedTables.Table(aoe = VectorOfSimilarVectors(aoe_grid), energy = e_rtft,
                blmean = bl_stats.mean, blslope = bl_stats.slope,
                t50 = t50_pre,
                qc_label = qc_labels
                )
end
export dsp_mwa_optimization_compressed
