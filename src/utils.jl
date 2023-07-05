7

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

    # fit window for basline extraction
    bl_mean::NTuple{2, Quantity{<:T}}

    # fit window for decay time extraction
    pz_fit::NTuple{2, Quantity{<:T}}

    # ADC threshold for t0 determination
    t0_threshold::T

    # rise and flat-top time grid scan ranges for trapezoidal filter
    e_grid_rt_trap::StepRangeLen{Quantity{<:T}, Base.TwicePrecision{Quantity{<:T}}, Base.TwicePrecision{Quantity{<:T}}, Int64}
    e_grid_ft_trap::StepRangeLen{Quantity{<:T}, Base.TwicePrecision{Quantity{<:T}}, Base.TwicePrecision{Quantity{<:T}}, Int64}

    # # rise and flat-top time grid scan ranges for ZAC filter
    e_grid_rt_zac::StepRangeLen{Quantity{<:T}, Base.TwicePrecision{Quantity{<:T}}, Base.TwicePrecision{Quantity{<:T}}, Int64}
    e_grid_ft_zac::StepRangeLen{Quantity{<:T}, Base.TwicePrecision{Quantity{<:T}}, Base.TwicePrecision{Quantity{<:T}}, Int64}

    # # rise and flat-top time grid scan ranges for CUSP filter
    e_grid_rt_cusp::StepRangeLen{Quantity{<:T}, Base.TwicePrecision{Quantity{<:T}}, Base.TwicePrecision{Quantity{<:T}}, Int64}
    e_grid_ft_cusp::StepRangeLen{Quantity{<:T}, Base.TwicePrecision{Quantity{<:T}}, Base.TwicePrecision{Quantity{<:T}}, Int64}
end
export DSPConfig

"""
    create_dsp_config(dsp_metadata::PropDicts.PropDict)

Create a `DSPConfig` from a `PropDict` of DSP metadata.

# Arguments
- `dsp_metadata::PropDicts.PropDict`: DSP metadata

# Returns
- `dsp_config::DSPConfig`: DSP configuration
"""
function create_dsp_config end
export create_dsp_config


function create_dsp_config(dsp_metadata::PropDicts.PropDict)
    return DSPConfig{Float64}(
        # pick-off time for ENC noise calculations
        dsp_metadata.enc_pickoff_trap*u"µs",
        dsp_metadata.enc_pickoff_zac*u"µs",
        dsp_metadata.enc_pickoff_cusp*u"µs",
        # fit window for basline extraction
        (dsp_metadata.bl_mean.min*u"µs", dsp_metadata.bl_mean.max*u"µs"),
        # fit window for decay time extraction
        (dsp_metadata.pz_fit.min*u"µs", dsp_metadata.pz_fit.max*u"µs"),
        # ADC threshold for t0 determination
        dsp_metadata.t0_threshold,
        # rise and flat-top time grid scan ranges for trapezoidal filter
        dsp_metadata.e_grid_trap.rt.start*u"µs":dsp_metadata.e_grid_trap.rt.step*u"µs":dsp_metadata.e_grid_trap.rt.stop*u"µs",
        dsp_metadata.e_grid_trap.ft.start*u"µs":dsp_metadata.e_grid_trap.ft.step*u"µs":dsp_metadata.e_grid_trap.ft.stop*u"µs",
        # rise and flat-top time grid scan ranges for ZAC filter
        dsp_metadata.e_grid_zac.rt.start*u"µs":dsp_metadata.e_grid_zac.rt.step*u"µs":dsp_metadata.e_grid_zac.rt.stop*u"µs",
        dsp_metadata.e_grid_zac.ft.start*u"µs":dsp_metadata.e_grid_zac.ft.step*u"µs":dsp_metadata.e_grid_zac.ft.stop*u"µs",
        # rise and flat-top time grid scan ranges for CUSP filter
        dsp_metadata.e_grid_cusp.rt.start*u"µs":dsp_metadata.e_grid_cusp.rt.step*u"µs":dsp_metadata.e_grid_cusp.rt.stop*u"µs",
        dsp_metadata.e_grid_cusp.ft.start*u"µs":dsp_metadata.e_grid_cusp.ft.step*u"µs":dsp_metadata.e_grid_cusp.ft.stop*u"µs")
end
