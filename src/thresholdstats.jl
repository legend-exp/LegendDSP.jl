# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    thresholdstats(signal::SamplesOrWaveform, min::RealQuantity=-Inf, max::RealQuantity=Inf)

Standard deviation of all samples of `signal` that lie between the lower bound `min` and the upper bound `max`.

The default bounds include all samples. The result is `NaN` if no sample lies between the bounds.
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
    thresholdstats_mad(signal::SamplesOrWaveform, min::RealQuantity=-Inf, max::RealQuantity=Inf)

Median absolute deviation (MAD) of all samples of `signal` that lie between the lower bound `min` and the upper bound `max`, scaled by 1.4826 to estimate the standard deviation of normally distributed data.

The default bounds include all samples. The result is zero if no sample lies between the bounds.
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
    @inbounds @fastmath for i in eachindex(Y_filt)
        Y_filt[i] = abs(Y_filt[i] - med)
    end
    # 1/Φ⁻¹(3/4) ≈ 1.4826, consistency constant for MAD under normality
    # (Rousseeuw & Croux, J. Amer. Statist. Assoc. 88, 1993)
    oftype(float(one(T)), 1.4826) * median!(Y_filt)
end