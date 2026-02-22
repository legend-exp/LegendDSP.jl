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
    thresholdstats_mad(signal, min=-Inf, max=Inf)

MAD-based noise estimate from mode: σ = 1.4826 × median(|Y - mode(Y)|).
Mode is approximated as rounded median (baseline level).
"""
function thresholdstats_mad(input::RadiationDetectorDSP.SamplesOrWaveform{T}, 
        min::RadiationDetectorDSP.RealQuantity = T(-Inf), 
        max::RadiationDetectorDSP.RealQuantity = T(Inf)) where T
    _, Y = RadiationDetectorDSP._get_axis_and_signal(input)
    _thresholdstats_mad_impl(Y, T(min), T(max))
end
export thresholdstats_mad

function _thresholdstats_mad_impl(Y::AbstractArray{T}, _min::T, _max::T) where {T <: RealQuantity}
    # Collect filtered samples
    Y_flt = T[]
    @inbounds for i in eachindex(Y)
        y = Y[i]
        if _min <= y <= _max
            push!(Y_flt, y)
        end
    end
    
    n = length(Y_flt)
    n == 0 && return zero(float(T))
    
    # Sort to find median
    sort!(Y_flt)
    med = if isodd(n)
        Y_flt[(n + 1) ÷ 2]
    else
        (Y_flt[n ÷ 2] + Y_flt[n ÷ 2 + 1]) / 2
    end
    
    # Mode approximation: rounded median (baseline level)
    mode_val = round(med; digits=0)
    
    # Compute absolute deviations from mode
    abs_dev = similar(Y_flt)
    @inbounds @fastmath @simd for i in eachindex(Y_flt)
        abs_dev[i] = abs(Y_flt[i] - mode_val)
    end
    
    # Sort to find median of absolute deviations
    sort!(abs_dev)
    mad = if isodd(n)
        abs_dev[(n + 1) ÷ 2]
    else
        (abs_dev[n ÷ 2] + abs_dev[n ÷ 2 + 1]) / 2
    end
    
    # MAD to σ conversion factor
    T(1.4826) * mad
end