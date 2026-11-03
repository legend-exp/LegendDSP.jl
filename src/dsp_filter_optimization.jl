# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    dsp_qc_flt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Unitful.Time{<:Real}, f_evaluate_qc::Union{Function, Missing})

Get QC DSP for filter parameter optimization for a given waveform set.
If f_evaluate_qc is missing, skips expensive ML classification and sets qc_label=-1.
Columns as for [`dsp_qc_flt_optimization_compressed`](@ref).
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

Columns: `energy` (default trapezoid), `blmean`, `blsigma`, `blslope`, `t50` (half of the
waveform maximum), `t50_pre` (half of the 10 µs rise, 10 µs flat-top trapezoid maximum, the `t50`
of `dsp_dplms_optimization`) and `qc_label`.
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

    # half of the current integral, as t50_pre in dsp_icpc_compressed and t50 in dsp_dplms_optimization
    e_101010 = extremestats.(TrapezoidalChargeFilter(10u"µs", 10u"µs").(wvfs)).max
    t50_pre = get_threshold(wvfs, e_101010 .* 0.5; mintot = config.kwargs_pars.tx_mintot)

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
                blmean = bl_stats.mean, blsigma = bl_stats.sigma, blslope = bl_stats.slope,
                t50 = t50, t50_pre = t50_pre,
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
    FltOptPrep

Baseline-subtracted, pole-zero corrected waveforms prepared once for repeated filter
evaluations. Created by [`dsp_flt_optimization_prepare`](@ref).
"""
struct FltOptPrep{W, T, B, E, P<:Unitful.Time{<:Real}}
    wvfs_pz::W
    t50::T
    bl_stats::B
    signal_estimator::E
    sampling_period::P

    FltOptPrep{W,T,B,E,P}(wvfs_pz, t50, bl_stats, signal_estimator, sampling_period) where {W,T,B,E,P} =
        new{W,T,B,E,P}(wvfs_pz, t50, bl_stats, signal_estimator, sampling_period)
end

FltOptPrep(wvfs_pz::W, t50::T, bl_stats::B, signal_estimator::E, sampling_period::P) where {W,T,B,E,P} =
    FltOptPrep{W,T,B,E,P}(wvfs_pz, t50, bl_stats, signal_estimator, sampling_period)
export FltOptPrep


"""
    dsp_flt_optimization_prepare(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Unitful.Time{<:Real})

Baseline-subtract and pole-zero correct `wvfs` once for repeated filter evaluations.
Returns a [`FltOptPrep`](@ref).
"""
function dsp_flt_optimization_prepare(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Unitful.Time{<:Real})
    wvfs = decode_data(wvfs)
    bl_stats = signalstats.(wvfs, leftendpoint(config.bl_window), rightendpoint(config.bl_window))
    wvfs_pz = InvCRFilter(τ).(shift_waveform.(wvfs, -bl_stats.mean))
    t50 = get_threshold(wvfs_pz, maximum.(wvfs_pz.signal) .* 0.5; mintot=config.kwargs_pars.tx_mintot)
    signal_estimator = SignalEstimator(PolynomialDNI(config.kwargs_pars.sig_interpolation_order, config.kwargs_pars.sig_interpolation_length))
    return FltOptPrep(wvfs_pz, t50, bl_stats, signal_estimator, step(first(wvfs).time))
end
export dsp_flt_optimization_prepare


"""
    dsp_energy_rtft(prep::FltOptPrep, config::DSPConfig, rt, ft, flt)

Energy of the prepared waveforms at rise time `rt` and flat-top time `ft`. `flt` is
`Val(:trap)`, `Val(:cusp)`, `Val(:zac)` or the corresponding `Symbol`.

The trapezoidal energy is read at `t50 + rt + ft/2`, the CUSP and ZAC energies at
`t50 + flt_length/2`. `rt` and `ft` are converted to floating point of the type of the
sampling period (trapezoidal) or the filter length (CUSP, ZAC). An energy is `NaN` if the
pickoff time lies outside the filtered time axis. Only reads `prep`, so it is safe to call
concurrently; each call allocates about two (trapezoidal) or five (CUSP, ZAC) times the
size of `prep.wvfs_pz`.
"""
function dsp_energy_rtft(prep::FltOptPrep, config::DSPConfig, rt::Unitful.Time{<:Real}, ft::Unitful.Time{<:Real}, ::Val{:trap})
    T = typeof(float(prep.sampling_period))
    uflt_rtft = TrapezoidalChargeFilter(convert(T, rt), convert(T, ft))
    return prep.signal_estimator.(uflt_rtft.(prep.wvfs_pz), prep.t50 .+ (rt + ft/2))
end

function dsp_energy_rtft(prep::FltOptPrep, config::DSPConfig, rt::Unitful.Time{<:Real}, ft::Unitful.Time{<:Real}, ::Val{:cusp})
    flt_length = config.flt_length_cusp
    scale = ustrip(NoUnits, flt_length / prep.sampling_period)
    T = typeof(float(flt_length))
    # τ = 10000000 µs disables the CR filter
    uflt_rtft = CUSPChargeFilter(convert(T, rt), convert(T, ft), convert(T, 10000000.0u"µs"), flt_length, scale)
    return prep.signal_estimator.(uflt_rtft.(prep.wvfs_pz), prep.t50 .+ flt_length/2)
end

function dsp_energy_rtft(prep::FltOptPrep, config::DSPConfig, rt::Unitful.Time{<:Real}, ft::Unitful.Time{<:Real}, ::Val{:zac})
    flt_length = config.flt_length_zac
    scale = ustrip(NoUnits, flt_length / prep.sampling_period)
    T = typeof(float(flt_length))
    uflt_rtft = ZACChargeFilter(convert(T, rt), convert(T, ft), convert(T, 10000000.0u"µs"), flt_length, scale)
    return prep.signal_estimator.(uflt_rtft.(prep.wvfs_pz), prep.t50 .+ flt_length/2)
end

dsp_energy_rtft(prep::FltOptPrep, config::DSPConfig, rt::Unitful.Time{<:Real}, ft::Unitful.Time{<:Real}, flt::Symbol) =
    dsp_energy_rtft(prep, config, rt, ft, Val(flt))
export dsp_energy_rtft

"""
    realized_flattop(ft, flt_length, sampling_period, flt)

Flat-top time that the filter `flt` (`Val(:trap)`, `Val(:cusp)`, `Val(:zac)` or the
corresponding `Symbol`) applies when built with flat-top time `ft` on waveforms sampled
with `sampling_period`. `flt_length` is the CUSP/ZAC filter length and is unused for `:trap`.

`RadiationDetectorDSP` rounds filter times to whole samples. The trapezoidal gap is
`round(ft / sampling_period)` samples. For CUSP and ZAC the flat top is lengthened by one
sample if the filter length minus the flat top, in samples, is odd. The result is a multiple
of `sampling_period`.
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
    dsp_trap_rt_optimization(prep::FltOptPrep, config::DSPConfig; ft::Quantity{T}=2.0u"µs")
    dsp_trap_rt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}; ft::Quantity{T}=2.0u"µs")

Get ENC noise values (read at `config.enc_pickoff_trap`) for each trap rise time in
`config.e_grid_rt_trap`. Returns a matrix of size `(n_rt, n_waveforms)`.
"""
function dsp_trap_rt_optimization(prep::FltOptPrep, config::DSPConfig; ft::Quantity{T}=2.0u"µs") where T<:Real
    enc_grid = zeros(Float64, length(config.e_grid_rt_trap), length(prep.wvfs_pz))
    for (r, rt) in enumerate(config.e_grid_rt_trap)
        uflt_rtft = TrapezoidalChargeFilter(rt, ft)
        enc_grid[r, :] = prep.signal_estimator.(uflt_rtft.(prep.wvfs_pz), config.enc_pickoff_trap)
    end
    return enc_grid
end

dsp_trap_rt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}; ft::Quantity{T}=2.0u"µs") where T<:Real =
    dsp_trap_rt_optimization(dsp_flt_optimization_prepare(wvfs, config, τ), config; ft)
export dsp_trap_rt_optimization


"""
    dsp_cusp_rt_optimization(prep::FltOptPrep, config::DSPConfig; ft::Quantity{T}=2.0u"µs")
    dsp_cusp_rt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}; ft::Quantity{T}=2.0u"µs")

Get ENC noise values (read at `config.enc_pickoff_cusp`) for each CUSP rise time in
`config.e_grid_rt_cusp`. Returns a matrix of size `(n_rt, n_waveforms)`.
"""
function dsp_cusp_rt_optimization(prep::FltOptPrep, config::DSPConfig; ft::Quantity{T}=2.0u"µs") where T<:Real
    flt_length = config.flt_length_cusp
    scale = ustrip(NoUnits, flt_length / prep.sampling_period)
    enc_grid = zeros(Float64, length(config.e_grid_rt_cusp), length(prep.wvfs_pz))
    for (r, rt) in enumerate(config.e_grid_rt_cusp)
        # τ = 10000000 µs disables the CR filter
        uflt_rtft = CUSPChargeFilter(rt, ft, 10000000.0u"µs", flt_length, scale)
        enc_grid[r, :] = prep.signal_estimator.(uflt_rtft.(prep.wvfs_pz), config.enc_pickoff_cusp)
    end
    return enc_grid
end

dsp_cusp_rt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}; ft::Quantity{T}=2.0u"µs") where T<:Real =
    dsp_cusp_rt_optimization(dsp_flt_optimization_prepare(wvfs, config, τ), config; ft)
export dsp_cusp_rt_optimization


"""
    dsp_zac_rt_optimization(prep::FltOptPrep, config::DSPConfig; ft::Quantity{T}=2.0u"µs")
    dsp_zac_rt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}; ft::Quantity{T}=2.0u"µs")

Get ENC noise values (read at `config.enc_pickoff_zac`) for each ZAC rise time in
`config.e_grid_rt_zac`. Returns a matrix of size `(n_rt, n_waveforms)`.
"""
function dsp_zac_rt_optimization(prep::FltOptPrep, config::DSPConfig; ft::Quantity{T}=2.0u"µs") where T<:Real
    flt_length = config.flt_length_zac
    scale = ustrip(NoUnits, flt_length / prep.sampling_period)
    enc_grid = zeros(Float64, length(config.e_grid_rt_zac), length(prep.wvfs_pz))
    for (r, rt) in enumerate(config.e_grid_rt_zac)
        uflt_rtft = ZACChargeFilter(rt, ft, 10000000.0u"µs", flt_length, scale)
        enc_grid[r, :] = prep.signal_estimator.(uflt_rtft.(prep.wvfs_pz), config.enc_pickoff_zac)
    end
    return enc_grid
end

dsp_zac_rt_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}; ft::Quantity{T}=2.0u"µs") where T<:Real =
    dsp_zac_rt_optimization(dsp_flt_optimization_prepare(wvfs, config, τ), config; ft)
export dsp_zac_rt_optimization


"""
    dsp_trap_ft_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{T}, rt::Quantity{T}) where T<:Real

Get trap energies at rise time `rt` for each flat-top time in `config.e_grid_ft_trap`.
Returns a matrix of size `(n_ft, n_waveforms)`.
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

Get CUSP energies at rise time `rt` for each flat-top time in `config.e_grid_ft_cusp`.
Returns a matrix of size `(n_ft, n_waveforms)`.
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

Get ZAC energies at rise time `rt` for each flat-top time in `config.e_grid_ft_zac`.
Returns a matrix of size `(n_ft, n_waveforms)`.
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

"""
    dplms_default_grid(n_rise::Integer; nm = 1.0, za = 1e3, ft = 1.0, ct = 0.0, zo = (0.0, 1e3), ff_factors = 1:4, ff = nothing)

DPLMS grid of flat-top widenings `ff = k·n_rise` samples for `k` in `ff_factors`, or of the
sample counts `ff` if given, for every
output-area weight in `zo` (a number or a collection), all at the weights `nm`, `za`, `ft`,
`ct` in units of the reference amplitude. Grid points are ordered by `zo`, then by `ff`. Only
weight ratios matter, and in these units any `za` above about 1 already enforces zero area;
`za = 0` switches the zero-area condition off. `zo = 0` leaves the area of the filter output
free, `zo = 1e3` makes it zero.
"""
function dplms_default_grid(n_rise::Integer; nm::Real = 1.0, za::Real = 1e3, ft::Real = 1.0, ct::Real = 0.0, zo = (0.0, 1e3), ff_factors = 1:4, ff = nothing)
    ffs = isnothing(ff) ? [round(Int, k * n_rise) for k in ff_factors] : ff
    [(nm = float(nm), za = float(za), ft = float(ft), ct = float(ct), zo = float(z), ff = Int(f)) for z in zo for f in ffs]
end
export dplms_default_grid

"""
    dsp_dplms_energy(wvfs::ArrayOfRDWaveforms, coeffs::AbstractVector, t50::AbstractVector, t_pickoff, current_window::ClosedInterval, signal_estimator::SignalEstimator)

Apply the DPLMS filter `coeffs` to baseline-subtracted, pole-zero-corrected waveforms and read
the energy at `t50 + t_pickoff` with `signal_estimator`. Returns `(e, emax)` with `emax` the
maximum of the filter output inside `current_window` shifted by `t_pickoff`, where the read-out
times of all events with `t50` inside `current_window` lie.
"""
function dsp_dplms_energy(wvfs::ArrayOfRDWaveforms, coeffs::AbstractVector, t50::AbstractVector, t_pickoff, current_window::ClosedInterval, signal_estimator::SignalEstimator)
    axes(wvfs) == axes(t50) || throw(DimensionMismatch("waveforms $(axes(wvfs)) and t50 times $(axes(t50)) differ"))
    wvfs = decode_data(wvfs)
    n_wvf = minimum(length, wvfs.signal)
    length(coeffs) <= n_wvf || throw(ArgumentError("$(length(coeffs)) DPLMS coefficients do not fit into waveforms of $n_wvf samples"))
    wvfs_flt = ConvolutionFilter(DirectConvolution(), coeffs).(wvfs)
    emax = get_wvf_maximum.(wvfs_flt, leftendpoint(current_window) + t_pickoff, rightendpoint(current_window) + t_pickoff)
    (e = signal_estimator.(wvfs_flt, t50 .+ t_pickoff), emax = [a.max for a in emax])
end
export dsp_dplms_energy

"""
    dsp_dplms_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity, drift::AbstractVector;
        n_pre = nothing, n_post = nothing, n_noise_windows = 8, ff_level = 0.999, n_bins = 5, plateau_length = 20, flat_tol = 1e-3, grid = dplms_default_grid)

Synthesize DPLMS energy filters on presummed waveforms for every point of `grid` and return the
energies they give. A grid point is a `NamedTuple` with the weights `nm`, `za`, `ft`, `ct`, `zo`
and the flat-top widening `ff` in samples.

Steps: baseline subtraction over `config.bl_window`; pole-zero correction with `InvCRFilter(τ)`;
`t50`, the time where half of the current integral is reached: 50 % of the 10/10 µs trapezoid
maximum, found with `get_threshold` (`config.kwargs_pars.tx_mintot`) as `t50_pre` in
`dsp_icpc_compressed`; alignment of every waveform to the sample nearest the median `t50`;
windows of `n_pre` samples before to `n_post` after that sample; noise matrix from the baseline
window; reference and one reference per drift-time bin, the bins being `n_bins`
quantiles of `drift`; flat-top matrix from the drift-bin references widened by `ff` samples.
All matrices are in units of the reference amplitude, the mean over the last `plateau_length`
window samples.

Unless given, `n_pre` and `n_post` equal `(n_bl - n_noise_windows) ÷ 2` with `n_bl` the number of
samples inside `config.bl_window`: the longest symmetric window that leaves `n_noise_windows`
sliding windows per event for the noise matrix. A window as long as the baseline leaves the
constant mode of the noise matrix undetermined, because the baseline mean is subtracted over the
same region, so only zero-area filters (`za > 0`) are meaningful for long windows.

Every grid point yields a filter applied to the unaligned pole-zero-corrected waveforms with
the energy read at `t50 + t_pickoff`, `t_pickoff = n_post·Δt`. The flat-top term
penalizes the output slope of every drift-bin reference shifted by up to `ff` samples, so the
output stays flat while any of these shapes moves by `ff` samples. `n_rise` is the number of
samples after the alignment point until the slowest drift-bin reference reaches `ff_level` of
the plateau. `grid` is either a vector of grid points or a function of `n_rise` returning one,
e.g. `grid = n -> dplms_default_grid(n; za = 0, zo = 0)` for filters without the zero-area condition.

`flat_top` and `shaping_time` describe each filter's output for the full-length drift-bin
references. The flat top is the interval around the read-out time where every filtered
reference stays within `flat_tol` of its value at the read-out time; the shaping time runs from
the first sample where the filtered superpulse exceeds `flat_tol` to the start of the flat top.

# Returns
- `e_grid`, `emax_grid`: energies and filter-output maxima, size `(length(grid), length(wvfs))`
- `coeffs`: impulse response per grid point; `t_pickoff`
- `noise_out`, `zero_area`: predicted relative output noise and coefficient sum per grid point
- `output_area`: area of the filter output for a unit step per grid point, in time units
- `flat_top`, `shaping_time`: one time per grid point
- `grid`, `ref`, `ref_bins` (one column per bin), `t50` (per event), `plateau`, `n_pre`, `n_post`, `n_rise`
- `superpulse`: the aligned superpulse in units of the reference amplitude
- `superpulse_bins`: the aligned superpulse of every drift-time bin, same units
- `response`: `superpulse` filtered with `coeffs`, one `RDWaveform` per grid point
- `noise`, `flattops`: the noise matrix and the flat-top matrix per `ff` (a `Dict`), in units of
  the reference amplitude
"""
function dsp_dplms_optimization(wvfs::ArrayOfRDWaveforms, config::DSPConfig, τ::Quantity{<:Real}, drift::AbstractVector;
        n_pre::Union{Integer,Nothing} = nothing, n_post::Union{Integer,Nothing} = nothing, n_noise_windows::Integer = 8,
        ff_level::Real = 0.999, n_bins::Integer = 5, plateau_length::Integer = 20,
        flat_tol::Real = 1e-3, grid = dplms_default_grid)
    axes(wvfs) == axes(drift) || throw(DimensionMismatch("waveforms $(axes(wvfs)) and drift proxies $(axes(drift)) differ"))
    all(isfinite, drift) || throw(ArgumentError("drift proxies must be finite, got $(count(!isfinite, drift)) non-finite values"))
    wvfs = decode_data(wvfs)
    t_rec = _shared_time_axis(wvfs)
    n_bl = _sample_index(t_rec, rightendpoint(config.bl_window)) - _sample_index(t_rec, leftendpoint(config.bl_window)) + 1
    n_half = (n_bl - n_noise_windows) ÷ 2
    n_pre = something(n_pre, n_half)
    n_post = something(n_post, n_half)
    n = n_pre + n_post + 1
    n <= n_bl || throw(ArgumentError("window of n_pre + n_post + 1 = $n samples is longer than the baseline of $n_bl samples"))
    1 <= plateau_length <= n_post || throw(ArgumentError("plateau_length $plateau_length must lie in 1:$n_post"))
    bl_window      = config.bl_window
    current_window = config.current_window

    bl_stats = signalstats.(wvfs, leftendpoint(bl_window), rightendpoint(bl_window))
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)
    wvfs = InvCRFilter(τ).(wvfs)
    Δt = step(first(wvfs.time))

    signal_estimator = SignalEstimator(PolynomialDNI(config.kwargs_pars.sig_interpolation_order, config.kwargs_pars.sig_interpolation_length))

    # half of the current integral, as t50_pre in dsp_icpc_compressed
    e_101010 = extremestats.(TrapezoidalChargeFilter(10u"µs", 10u"µs").(wvfs)).max
    t50 = get_threshold(wvfs, e_101010 .* 0.5; mintot = config.kwargs_pars.tx_mintot)
    n_out = count(∉(current_window), t50)
    n_out == 0 || throw(ArgumentError("t50 lies outside the current window $current_window for $n_out of $(length(t50)) waveforms"))

    # aligned windows in units of the reference amplitude
    t = first(wvfs.time)
    t_target = t[_sample_index(t, median(t50))]
    wvfs_aligned = align_waveforms(wvfs, t50, t_target)
    windows = aligned_windows(wvfs_aligned, t_target, n_pre, n_post)
    ref, plateau = dplms_reference(windows, (n - plateau_length + 1):n)
    noise = dplms_noise_matrix(wvfs, bl_window, n) ./ plateau^2
    ranks = invperm(sortperm(drift))
    bins = cld.(ranks .* n_bins, length(drift))
    binned = dplms_binned_reference(windows, bins, plateau)
    ref_bins = binned.refs

    # samples from the alignment point until the slowest bin reaches ff_level
    n_rise = maximum(axes(ref_bins, 2)) do b
        i = findfirst(>=(ff_level), view(ref_bins, (n_pre + 1):n, b))
        isnothing(i) && throw(ArgumentError("drift-time bin $b never reaches $ff_level of the plateau after the alignment point"))
        max(i - 1, 1)
    end
    grid = grid isa AbstractVector ? grid : grid(n_rise)
    ff_max = maximum(w -> w.ff, grid)
    refs_wide = dplms_binned_reference(aligned_windows(wvfs_aligned, t_target, n_pre + ff_max, n_post + ff_max), bins, plateau).refs
    flattops = Dict(ff => dplms_flattop_matrix(refs_wide[(1 + ff_max - ff):(end - ff_max + ff), :], ff) for ff in unique(w.ff for w in grid))
    superpulse = RDWaveform(t, mean(wvfs_aligned.signal) ./ plateau)
    superpulse_bins = [RDWaveform(t, mean(wvfs_aligned.signal[bins .== b]) ./ plateau) for b in axes(ref_bins, 2)]

    t_pickoff = n_post * Δt
    e_grid    = Matrix{Float64}(undef, length(grid), length(wvfs))
    emax_grid = similar(e_grid)
    coeffs    = Vector{Vector{Float64}}(undef, length(grid))
    noise_out = Vector{Float64}(undef, length(grid))
    zero_area = similar(noise_out)
    output_area = Vector{typeof(float(Δt))}(undef, length(grid))
    for (g, w) in enumerate(grid)
        s = dplms_synthesis(noise, ref, flattops[w.ff], binned; nm = w.nm, za = w.za, ft = w.ft, ct = w.ct, zo = w.zo)
        energies = dsp_dplms_energy(wvfs, s.coeffs, t50, t_pickoff, current_window, signal_estimator)
        e_grid[g, :]    = energies.e
        emax_grid[g, :] = energies.emax
        coeffs[g]    = s.coeffs
        noise_out[g] = s.noise_out
        zero_area[g] = s.zero_area
        output_area[g] = s.output_area * Δt
    end
    response = [ConvolutionFilter(DirectConvolution(), c)(superpulse) for c in coeffs]

    shapes = map(coeffs, response) do c, resp
        outs = [ConvolutionFilter(DirectConvolution(), c)(sp).signal for sp in superpulse_bins]
        i0 = _sample_index(resp.time, t_target + t_pickoff)
        flat(i) = all(o -> abs(o[i] - o[i0]) <= flat_tol, outs)
        i_lo, i_hi = i0, i0
        while i_lo > firstindex(resp.signal) && flat(i_lo - 1); i_lo -= 1; end
        while i_hi < lastindex(resp.signal) && flat(i_hi + 1); i_hi += 1; end
        i_start = findfirst(>(flat_tol), resp.signal)
        (flat_top = (i_hi - i_lo) * Δt, shaping_time = (i_lo - something(i_start, i_lo)) * Δt)
    end
    flat_top = [x.flat_top for x in shapes]
    shaping_time = [x.shaping_time for x in shapes]

    (; e_grid, emax_grid, coeffs, t_pickoff, flat_top, shaping_time, noise_out, zero_area, output_area, grid, ref, ref_bins, t50, plateau, n_pre, n_post, n_rise, superpulse, superpulse_bins, response, noise, flattops)
end
export dsp_dplms_optimization
