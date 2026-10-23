# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    DSPConfig{T <: Real}

Configuration parameters for DSP algorithms.

# Fields
- `enc_pickoff::Quantity{<:T}`: pick-off time for ENC noise calculations
- `bl_window::ClosedInterval{Quantity{<:T}}`: fit window for basline extraction
- `tail_window::ClosedInterval{Quantity{<:T}}`: fit window for decay time extraction
- `inTraceCut_std_threshold::T`: in-trace pile-up rejector threshold in standard deviations
- `t0_threshold::T`: ADC threshold for t0 determination
- `e_grid_rt_trap::StepRangeLen{Quantity{<:T}}`: rise time grid scan range for trapezoidal filter
- `e_grid_ft_trap::StepRangeLen{Quantity{<:T}}`: flat-top time grid scan range for trapezoidal filter
- `e_grid_rt_zac::StepRangeLen{Quantity{<:T}}`: rise time grid scan range for ZAC filter
- `e_grid_ft_zac::StepRangeLen{Quantity{<:T}}`: flat-top time grid scan range for ZAC filter
- `e_grid_rt_cusp::StepRangeLen{Quantity{<:T}}`: rise time grid scan range for CUSP filter
- `e_grid_ft_cusp::StepRangeLen{Quantity{<:T}}`: flat-top time grid scan range for CUSP filter
- `a_grid_wl_sg::StepRangeLen{Quantity{<:T}}`: window length grid scan range for SG filter
- `a_grid_wl_mwa::StepRangeLen{Quantity{<:T}}`: window length grid scan range for MWA filter
- `a_mwa_windows::Vector{Quantity{<:T}}`: fixed MWA window lengths stored as `a_mwa_<ns>` columns
- `a_mwa_passes::Int`: number of alternating moving average passes of the MWA current filter
- `a_sg_windows::Vector{Quantity{<:T}}`: fixed SG window lengths stored as `a_sg_<ns>` columns
- `lq_windows::Vector{Vector{Quantity{<:T}}}`: `[early, gap, late]` integration windows of the LQ parameters; one column `lq80_<early><gap><late>` (and `lq90_…`, `lqamax_…`) per window

# Examples
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


function DSPConfig(pd::PropDict)
    _create_dsp_config(pd)
end


