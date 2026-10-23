# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    DSPConfig{T <: Real}

Settings of the HPGe DSP chains and filter-parameter scans.

All time-valued fields are `Unitful` quantities. Create an instance from DSP metadata with `DSPConfig(pd::PropDict)`.

# Fields
- `enc_pickoff_trap`, `enc_pickoff_zac`, `enc_pickoff_cusp`: time at which the trapezoidal, ZAC and CUSP filtered waveforms are evaluated in the rise-time scans [`dsp_trap_rt_optimization`](@ref), [`dsp_zac_rt_optimization`](@ref) and [`dsp_cusp_rt_optimization`](@ref)
- `flt_length_cusp`, `flt_length_zac`: total length of the CUSP and ZAC filters
- `t0_threshold`: threshold in units of the filtered signal for the t0 determination of [`get_t0`](@ref)
- `inTraceCut_std_threshold`: threshold of the in-trace pile-up search in standard deviations of the baseline
- `sg_flt_degree`: polynomial degree of the Savitzky-Golay filters that determine the current
- `bl_window`: time window for the baseline fit
- `tail_window`: time window for the decay-time fit and the tail statistics after pole-zero correction
- `current_window`: time window in which the maximum of the current is searched
- `qdrift_int_length`: `StepRangeLen` whose first and last elements are the two integration times of [`get_qdrift`](@ref)
- `lq_windows`: vector of `[early, gap, late]` window lengths of the LQ parameters, see [`get_lq`](@ref)
- `e_grid_rt_trap`, `e_grid_ft_trap`: rise-time and flat-top-time scan ranges of the trapezoidal filter
- `e_grid_rt_zac`, `e_grid_ft_zac`: rise-time and flat-top-time scan ranges of the ZAC filter
- `e_grid_rt_cusp`, `e_grid_ft_cusp`: rise-time and flat-top-time scan ranges of the CUSP filter
- `a_grid_wl_sg`: window-length scan range of the Savitzky-Golay current filter
- `a_grid_wl_mwa`: window-length scan range of the moving-window-average current filter
- `a_mwa_windows`: fixed window lengths of the moving-window-average current filter; each gives a column `a_mwa_<ns>`
- `a_mwa_passes`: number of alternating passes of the moving-window-average current filter, see [`MovingWindowMultiFilter`](@ref)
- `a_sg_windows`: fixed window lengths of the Savitzky-Golay current filter; each gives a column `a_sg_<ns>`
- `default_flt_param`: `PropDict` with the default filter parameters: `rt` and `ft` for `trap`, `cusp` and `zac`, and the window length for `sg` and `mwa`
- `kwargs_pars`: `PropDict` with further parameters: `fc_bit_depth`, `t0_flt_pars`, `t0_mintot`, `tx_mintot`, `intrace_mintot`, `int_interpolation_order`, `int_interpolation_length`, `sig_interpolation_order` and `sig_interpolation_length`
- `auxpz1_window`, `auxpz2_window`: time windows of the auxiliary pole-zero statistics `auxpz1_*` and `auxpz2_*` of [`dsp_icpc_compressed`](@ref)

# Examples
`LegendDataManagement` is not a dependency of `LegendDSP` and is only used here to read the metadata.
```julia
using LegendDSP
using LegendDataManagement

l200 = LegendData(:l200)
filekey = start_filekey(l200, (:p03, :r000, :cal))
dsp_config = DSPConfig(dataprod_config(l200).dsp(filekey).default)
```
"""
struct DSPConfig{T <: Real}
    # pick-off times for ENC noise calculations
    enc_pickoff_trap::Quantity{<:T}
    enc_pickoff_zac::Quantity{<:T}
    enc_pickoff_cusp::Quantity{<:T}

    # filter lengths for CUSP and ZAC filters
    flt_length_cusp::Quantity{<:T}
    flt_length_zac::Quantity{<:T}

    # ADC threshold for t0 determination
    t0_threshold::T

    # in-trace pile-up rejector threshold in sigmas
    inTraceCut_std_threshold::T

    # Savitzky Golay polynominal order for current extraction
    sg_flt_degree::Int

    # fit window for basline extraction
    bl_window::ClosedInterval{<:Quantity{<:T}}

    # fit window for decay time extraction
    tail_window::ClosedInterval{<:Quantity{<:T}}

    # fit window for current extraction
    current_window::ClosedInterval{<:Quantity{<:T}}

    # Integration length for QDrift extraction
    qdrift_int_length::StepRangeLen{<:Quantity{<:T}}

    # (early, gap, late) integration windows of the LQ parameters
    lq_windows::Vector{<:AbstractVector{<:Quantity{<:T}}}

    # rise and flat-top time grid scan ranges for trapezoidal filter
    e_grid_rt_trap::StepRangeLen{Quantity{<:T}, Base.TwicePrecision{Quantity{<:T}}, Base.TwicePrecision{Quantity{<:T}}, Int64}
    e_grid_ft_trap::StepRangeLen{Quantity{<:T}, Base.TwicePrecision{Quantity{<:T}}, Base.TwicePrecision{Quantity{<:T}}, Int64}

    # rise and flat-top time grid scan ranges for ZAC filter
    e_grid_rt_zac::StepRangeLen{Quantity{<:T}, Base.TwicePrecision{Quantity{<:T}}, Base.TwicePrecision{Quantity{<:T}}, Int64}
    e_grid_ft_zac::StepRangeLen{Quantity{<:T}, Base.TwicePrecision{Quantity{<:T}}, Base.TwicePrecision{Quantity{<:T}}, Int64}

    # rise and flat-top time grid scan ranges for CUSP filter
    e_grid_rt_cusp::StepRangeLen{Quantity{<:T}, Base.TwicePrecision{Quantity{<:T}}, Base.TwicePrecision{Quantity{<:T}}, Int64}
    e_grid_ft_cusp::StepRangeLen{Quantity{<:T}, Base.TwicePrecision{Quantity{<:T}}, Base.TwicePrecision{Quantity{<:T}}, Int64}
    
    # window length grid scan range for SG filter in current determination
    a_grid_wl_sg::StepRangeLen{Quantity{<:T}, Base.TwicePrecision{Quantity{<:T}}, Base.TwicePrecision{Quantity{<:T}}, Int64}

    # window length grid scan range for MWA filter in current determination
    a_grid_wl_mwa::StepRangeLen{Quantity{<:T}, Base.TwicePrecision{Quantity{<:T}}, Base.TwicePrecision{Quantity{<:T}}, Int64}

    # fixed MWA window lengths stored as `a_mwa_<ns>` columns and number of alternating moving average passes
    a_mwa_windows::Vector{<:Quantity{<:T}}
    a_mwa_passes::Int

    # fixed SG window lengths stored as `a_sg_<ns>` columns
    a_sg_windows::Vector{<:Quantity{<:T}}

    # default filter parameter
    default_flt_param::PropDict

    # additional pars
    kwargs_pars::PropDict

    # auxiliary pz-corrected waveform windows
    auxpz1_window::ClosedInterval{<:Quantity{<:T}}
    auxpz2_window::ClosedInterval{<:Quantity{<:T}}
end
export DSPConfig


"""
    DSPConfig(pd::PropDict)

Create a `DSPConfig{Float64}` from DSP metadata.

`pd` needs one entry per field of [`DSPConfig`](@ref), with these differences:

- `bl_window`, `tail_window`, `current_window`, `auxpz1_window` and `auxpz2_window` are given as `(min, max)` entries;
- `e_grid_trap`, `e_grid_zac` and `e_grid_cusp` each contain `rt` and `ft` entries with `start`, `step` and `stop`; `a_grid_wl_sg` and `a_grid_wl_mwa` have the same three entries. `a_grid_wl_mwa` is optional and defaults to `a_grid_wl_sg`;
- `qdrift_int_length` holds a start and a stop time; the resulting range has a step of 0.1 in the unit of the start time;
- `default_flt_param` is read from the entry `flt_defaults`;
- every element of `lq_windows` must have exactly three entries, otherwise an `ArgumentError` is thrown.
"""
function DSPConfig(pd::PropDict)
    _create_dsp_config(pd)
end


