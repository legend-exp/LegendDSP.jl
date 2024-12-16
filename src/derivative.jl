# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    struct DerivativeFilter <: AbstractRadIIRFilter

Working example:
```julia
using LegendDSP
using Unitful

signal = rand(100)
t = range(0u"ms", 20u"ms", 100)

wf = RDWaveform(t, signal)
flt = DerivativeFilter()
wf_new = flt(wf)
```
Constructors:

* ```$(FUNCTIONNAME)(; fields...)```

Fields:

$(TYPEDFIELDS)
"""
Base.@kwdef struct DerivativeFilter{T <: RadiationDetectorDSP.RealQuantity} <: AbstractRadIIRFilter
    "Filter gain"
    gain::T = 1
end
export DerivativeFilter

struct DerivativeFilterInstance{T <: RadiationDetectorDSP.RealQuantity, TG <: RadiationDetectorDSP.RealQuantity} <: AbstractRadSigFilterInstance{LinearFiltering}
    gain::TG
    n_input::Int
end

function RadiationDetectorDSP.fltinstance(flt::DerivativeFilter{TG}, si::SamplingInfo{T}) where {T, TG}
    DerivativeFilterInstance{T, TG}(flt.gain, _smpllen(si))
end

RadiationDetectorDSP.flt_output_smpltype(::DerivativeFilterInstance{T, TG})  where {T,TG} = promote_type(typeof(zero(T) * unit(TG)), TG)
RadiationDetectorDSP.flt_input_smpltype(::DerivativeFilterInstance{T}) where {T} = T
RadiationDetectorDSP.flt_output_length(fi::DerivativeFilterInstance) = fi.n_input
RadiationDetectorDSP.flt_input_length(fi::DerivativeFilterInstance) = fi.n_input
RadiationDetectorDSP.flt_output_time_axis(::DerivativeFilterInstance, time::AbstractVector) = time

function rdfilt!(y::AbstractVector, fi::DerivativeFilterInstance, x::AbstractVector)

    @assert firstindex(y) == firstindex(x)
    @assert lastindex(y) == lastindex(x)
    @inbounds @simd for i in eachindex(y)
        y[i] = fi.gain * (x[max(i,begin+1)] - x[max(i-1,begin)])
    end
    y
end