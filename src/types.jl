# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

const MaybeWithUnits{T<:Number} = Union{T,Quantity{<:T}}
const RealQuantity = MaybeWithUnits{<:Real}


"""
    DSPConfig{T <: Real}

Configuration parameters for DSP algorithms.

# Fields
- `enc_pickoff::Quantity{<:T}`: pick-off time for ENC noise calculations
- `bl_mean::NTuple{2, Quantity{<:T}}`: fit window for basline extraction
- `pz_fit::NTuple{2, Quantity{<:T}}`: fit window for decay time extraction
- `t0_threshold::T`: ADC threshold for t0 determination
- `e_grid_rt_trap::StepRangeLen{Quantity{<:T}}`: rise time grid scan range for trapezoidal filter
- `e_grid_ft_trap::StepRangeLen{Quantity{<:T}}`: flat-top time grid scan range for trapezoidal filter
- `e_grid_rt_zac::StepRangeLen{Quantity{<:T}}`: rise time grid scan range for ZAC filter
- `e_grid_ft_zac::StepRangeLen{Quantity{<:T}}`: flat-top time grid scan range for ZAC filter
- `e_grid_rt_cusp::StepRangeLen{Quantity{<:T}}`: rise time grid scan range for CUSP filter
- `e_grid_ft_cusp::StepRangeLen{Quantity{<:T}}`: flat-top time grid scan range for CUSP filter
- `e_grid_rt_sg::StepRangeLen{Quantity{<:T}}`: window length grid scan range for SG filter

# Examples
```julia
using LegendDSP
using Unitful

dsp_config = DSPConfig{Float64}(32.0u"µs", 
(0.0u"µs", 39.0u"µs"), 
(80.0u"µs", 110.0u"µs"), 
5.0, 
7.0u"µs":0.5u"µs":12.0u"µs", 1.0u"µs":0.2u"µs":4.0u"µs",
7.0u"µs":0.5u"µs":12.0u"µs", 1.0u"µs":0.2u"µs":4.0u"µs", 
7.0u"µs":0.5u"µs":12.0u"µs", 1.0u"µs":0.2u"µs":4.0u"µs")
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

    # in-trace pile-up rejector threshold in standard deviations
    inTraceCut_std_threshold::T

    # fit window for basline extraction
    bl_mean::NTuple{2, Quantity{<:T}}

    # fit window for decay time extraction
    pz_fit::NTuple{2, Quantity{<:T}}

    # ADC threshold for t0 determination
    t0_threshold::T

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
end
export DSPConfig



