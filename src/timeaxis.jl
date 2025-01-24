# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    struct TimeAxisFilter <: AbstractRadIIRFilter

Updates the time range step to `period` and shifts the time axis by `offset`.

Working example:
```julia
using LegendDSP
using RadiationDetectorSignals
using Unitful

signal = rand(100)
t = range(0u"ms", 20u"ms", 100)

wf = RDWaveform(t, signal)
flt = TimeAxisFilter(4u"ns")
wf_new = flt(wf)
```
Constructors:

* ```$(FUNCTIONNAME)(; fields...)```

Fields:

$(TYPEDFIELDS)
"""
Base.@kwdef struct TimeAxisFilter{T <: RadiationDetectorDSP.RealQuantity} <: AbstractRadIIRFilter
    "Filter period"
    period::T
    "Filter offset"
    offset::T = false * period
    function TimeAxisFilter(period::T, offset::T=false * period) where T
        new{T}(period, offset)
    end
end
export TimeAxisFilter

struct TimeAxisFilterInstance{T <: RadiationDetectorDSP.RealQuantity, TT <: RadiationDetectorDSP.RealQuantity} <: AbstractRadSigFilterInstance{LinearFiltering}
    period::TT
    offset::TT
    n_input::Int
end

function RadiationDetectorDSP.fltinstance(flt::TimeAxisFilter{TT}, si::SamplingInfo{T}) where {T <: RadiationDetectorDSP.RealQuantity, TT <: RadiationDetectorDSP.RealQuantity}
    TimeAxisFilterInstance{T,TT}(flt.period, flt.offset, _smpllen(si))
end

RadiationDetectorDSP.flt_output_smpltype(fi::TimeAxisFilterInstance{T})  where {T} = flt_input_smpltype(fi)
RadiationDetectorDSP.flt_input_smpltype(::TimeAxisFilterInstance{T}) where {T} = T
RadiationDetectorDSP.flt_output_length(fi::TimeAxisFilterInstance) = fi.n_input
RadiationDetectorDSP.flt_input_length(fi::TimeAxisFilterInstance) = fi.n_input
RadiationDetectorDSP.flt_output_time_axis(fi::TimeAxisFilterInstance, time::AbstractVector) = throw(ArgumentError("Non-range time axis not supported"))
RadiationDetectorDSP.flt_output_time_axis(fi::TimeAxisFilterInstance, time::AbstractRange) = range(first(time) + fi.offset, step=fi.period, length=length(time))

rdfilt!(y::AbstractVector, fi::TimeAxisFilterInstance, x::AbstractVector) = y .= x
