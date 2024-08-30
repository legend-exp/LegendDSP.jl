import RadiationDetectorDSP: AbstractConvFilterInstance,
    fltinstance, _smpllen, smplinfo, SamplingInfo, 
    AbstractRadSigFilterInstance, LinearFiltering, _filterlen, _floattype

"""
    struct MovingWindowMulti <: AbstractRadFIRFilter

apply left and right moving average windows to signal.

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
    invl = 1/l
    @inbounds @simd for i in 2:l
        y[i] = fma(invl, x[i] - x₁, y[i-1])     # yₙ = yₙ₋₁ + (xₙ - x₁)/l
    end
    L = length(y)
    @inbounds @simd for i in l+1:L
        y[i] = fma(invl, x[i] - x[i - l], y[i-1]) # yₙ = yₙ₋₁ + (xₙ - xₙ₋ₗ)/l
    end
    y
end

function rdfilt!(y::AbstractVector, fi::MovingWindowMultiFilterInstance{T}, x::AbstractVector
    ) where T

    _fi = MovingWindowFilterInstance{T}(fi.length, fi.n_input)
    _y = similar(y)
    yr = @view y[end:-1:begin]          # reversed view of y
    _yr = @view _y[end:-1:begin]        # reversed view of _y
    rdfilt!(_y, _fi, x)
    rdfilt!(yr, _fi, _yr)
    copyto!(_y, y)
    rdfilt!(y, _fi, _y)    
end