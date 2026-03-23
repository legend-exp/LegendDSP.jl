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


"""
    thresholdstats_mad(signal::AbstractSamples, min::Real, max::Real)
    thresholdstats_mad(signal::RDWaveform, min::RealQuantity, max::RealQuantity)

Get a robust threshold estimate using the Median Absolute Deviation (MAD)
of all samples in a `signal` that are within the lower bound `min` and
the upper bound `max`. The MAD is scaled by 1.4826 to be consistent with
the standard deviation for normally distributed data.
"""
function thresholdstats_mad end
export thresholdstats_mad

function thresholdstats_mad(input::RadiationDetectorDSP.SamplesOrWaveform{T}, min::RadiationDetectorDSP.RealQuantity = -Inf * unit(T), max::RadiationDetectorDSP.RealQuantity = Inf * unit(T)) where T
    _, Y = RadiationDetectorDSP._get_axis_and_signal(input)
    _thresholdstats_mad_impl(Y, T(min), T(max))
end

function _thresholdstats_mad_impl(Y::AbstractArray{T}, _min::T, _max::T) where {T <: RealQuantity}
    Y_filt = [y for y in Y if _min <= y <= _max]
    isempty(Y_filt) && return zero(float(T))
    med = median(Y_filt)
    T(1.4826) * median(abs.(Y_filt .- med))
end