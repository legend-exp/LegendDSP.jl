# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    thresholdstats(signal::AbstractSamples, min::Real, min::Real)
    thresholdstats(signal::RDWaveform, min::RealQuantity, max::RealQuantity)

Get the standard deviation of all samples in a `signal` that are within 
the lower bound `min` and the upper bound `max`.
If no values for `min` or `max` are passed, the respective bound is ignored.
"""
function thresholdstats end
export thresholdstats

function thresholdstats(input::RadiationDetectorDSP.SamplesOrWaveform{T}, min::RadiationDetectorDSP.RealQuantity = -Inf * unit(T), max::RadiationDetectorDSP.RealQuantity = Inf * unit(T)) where T
    _, Y = RadiationDetectorDSP._get_axis_and_signal(input)
    _thresholdstats_impl(Y, T(min), T(max))
end

function _thresholdstats_impl(Y::AbstractArray{T}, _min::T, _max::T) where {T <: RealQuantity}
        
    zy = zero(eltype(Y))

    sum_Y::float(typeof(zy)) = zy
    sum_Y_sqr::float(typeof(zy * zy)) = zy * zy
    n::Int = 0

    @inbounds @fastmath @simd for i in eachindex(Y)
        y = Y[i]
        _include = _min <= y <= _max
        y = _include * y
        sum_Y = y + sum_Y
        sum_Y_sqr = fma(y, y, sum_Y_sqr)
        n += _include
    end
    
    inv_n = inv(n)
    mean_Y = sum_Y * inv_n
    var_Y = max(sum_Y_sqr * inv_n - mean_Y * mean_Y, zero(sum_Y_sqr))
    sqrt(var_Y)

end