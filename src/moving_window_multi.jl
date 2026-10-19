import RadiationDetectorDSP: AbstractConvFilterInstance,
    fltinstance, _smpllen, smplinfo, SamplingInfo, 
    AbstractRadSigFilterInstance, LinearFiltering, _filterlen, _floattype

"""
    struct MovingWindowMultiFilter <: AbstractRadFIRFilter

apply left and right moving average windows to signal.
Working example:
```julia
using RadiationDetectorSignals
using Unitful

signal = rand(100)
t = range(0u"ms", 20u"ms", 100)

wf = RDWaveform(t, signal)
flt = MovingWindowMultiFilter(1u"ms")
wf_new = flt(wf)
```
Constructors:

* ```$(FUNCTIONNAME)(; fields...)```

Fields:

$(TYPEDFIELDS)
"""
Base.@kwdef struct MovingWindowMultiFilter{T <: RadiationDetectorDSP.RealQuantity} <: AbstractRadFIRFilter
    "size of the moving average window"
    length::T
end

export MovingWindowMultiFilter

struct MovingWindowMultiFilterInstance{T} <: AbstractRadSigFilterInstance{LinearFiltering}
    length::Int
    n_input::Int
end

"""
    struct MovingWindowFilter <: AbstractRadFIRFilter

apply left moving average window to signal. The exact computations
are:\\
        yₙ = yₙ₋₁ + (xₙ - x₁)/l,    for n ∈ {1, …, l}\\
        yₙ = yₙ₋₁ + (xₙ - xₙ₋ₗ)/l,    for n ∈ {l+1, …, L}\\

    
Constructors:

* ```$(FUNCTIONNAME)(; fields...)```

Fields:

$(TYPEDFIELDS)
"""
Base.@kwdef struct MovingWindowFilter{T <: RadiationDetectorDSP.RealQuantity} <: AbstractRadFIRFilter
    "length of the moving average window"
    length::T
end

export MovingWindowFilter

struct MovingWindowFilterInstance{T} <: AbstractRadSigFilterInstance{LinearFiltering}
    length::Int
    n_input::Int
end

function RadiationDetectorDSP.fltinstance(flt::MovingWindowFilter, 
    si::SamplingInfo{T}) where {T}

    length = round(Int, uconvert(NoUnits, flt.length / step(si.axis)))
    MovingWindowFilterInstance{T}(length, length + _smpllen(si) - 1)
end

function RadiationDetectorDSP.fltinstance(flt::MovingWindowMultiFilter, 
    si::SamplingInfo{T}) where {T}

    length = round(Int, uconvert(NoUnits, flt.length / step(si.axis)))
    MovingWindowMultiFilterInstance{T}(length, length + _smpllen(si) - 1)
end

RadiationDetectorDSP._filterlen(fi::MovingWindowFilterInstance) = fi.length
RadiationDetectorDSP._filterlen(fi::MovingWindowMultiFilterInstance) = fi.length

RadiationDetectorDSP.flt_output_smpltype(::MovingWindowFilterInstance{T})  where T = _floattype(T)
RadiationDetectorDSP.flt_input_smpltype(::MovingWindowFilterInstance{T}) where {T} = T
RadiationDetectorDSP.flt_output_length(fi::MovingWindowFilterInstance) = RadiationDetectorDSP.flt_input_length(fi) - _filterlen(fi) + 1
RadiationDetectorDSP.flt_input_length(fi::MovingWindowFilterInstance) = fi.n_input
RadiationDetectorDSP.flt_output_time_axis(::MovingWindowFilterInstance, time::AbstractVector) = time

RadiationDetectorDSP.flt_output_smpltype(::MovingWindowMultiFilterInstance{T}) where {T} = _floattype(T)
RadiationDetectorDSP.flt_input_smpltype(::MovingWindowMultiFilterInstance{T}) where {T} = T
RadiationDetectorDSP.flt_output_length(fi::MovingWindowMultiFilterInstance) = RadiationDetectorDSP.flt_input_length(fi) - _filterlen(fi) + 1
RadiationDetectorDSP.flt_input_length(fi::MovingWindowMultiFilterInstance) = fi.n_input
RadiationDetectorDSP.flt_output_time_axis(::MovingWindowMultiFilterInstance, time::AbstractVector) = time

function rdfilt!(y::AbstractVector, fi::MovingWindowFilterInstance{T}, x::AbstractVector
    ) where T

    @assert firstindex(y) == firstindex(x)
    @assert lastindex(y) == lastindex(x)
    x₁ = x[begin]
    y[begin] = x₁
    l = fi.length
    invl = one(eltype(y)) / l
    @inbounds @simd for i in 2:l
        y[i] = fma(invl, x[i] - x₁, y[i-1])     # yₙ = yₙ₋₁ + (xₙ - x₁)/l
    end
    L = length(y)
    @inbounds @simd for i in l+1:L
        y[i] = fma(invl, x[i] - x[i - l], y[i-1]) # yₙ = yₙ₋₁ + (xₙ - xₙ₋ₗ)/l
    end
    y
end

function rdfilt!(y::AbstractVector, fi::MovingWindowMultiFilterInstance, x::AbstractVector)
    axes(y) == axes(x) || throw(DimensionMismatch("input and output must have the same axes, got $(axes(x)) and $(axes(y))"))
    l = fi.length
    l >= 1 || throw(ArgumentError("moving window length must be at least one sample, got $l"))
    # moving window forward, backward and forward again; the second and third pass work in place on y
    _moving_window_pass!(y, x, eachindex(y), l, Vector{eltype(x)}(undef, l))
    buf = Vector{eltype(y)}(undef, l)
    _moving_window_pass!(y, y, reverse(eachindex(y)), l, buf)
    _moving_window_pass!(y, y, eachindex(y), l, buf)
end

# One pass of the MovingWindowFilter recursion over the indices `idxs` (in this order):
#     y[idxs[1]] = x[idxs[1]],  y[idxs[k]] = y[idxs[k-1]] + (x[idxs[k]] - x[idxs[max(k-l, 1)]]) / l
# The ring buffer `buf` (length `l`) holds the last `l` input values, so `y` may be the same array as `x`.
function _moving_window_pass!(y::AbstractVector, x::AbstractVector, idxs, l::Int, buf::AbstractVector)
    invl = one(eltype(y)) / l
    i₁ = first(idxs)
    x₁ = x[i₁]
    acc = convert(eltype(y), x₁)
    y[i₁] = acc
    # input values before the first index count as x₁
    fill!(buf, x₁)
    j = firstindex(buf)
    for i in Iterators.drop(idxs, 1)
        xᵢ = x[i]
        acc = convert(eltype(y), fma(invl, xᵢ - buf[j], acc))
        y[i] = acc
        buf[j] = xᵢ
        j = j == lastindex(buf) ? firstindex(buf) : j + 1
    end
    y
end