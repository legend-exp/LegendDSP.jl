# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    saturation(signal::SamplesOrWaveform, low::RealQuantity, high::RealQuantity)
    saturation(signal::SamplesOrWaveform, start::RealQuantity, stop::RealQuantity, low::RealQuantity, high::RealQuantity)

Count the samples of `signal` that are equal to the lower saturation value `low` or the upper saturation value `high`.

The first method uses the whole signal, the second the interval between `start` and `stop`, including the samples nearest to both.

# Returns
- named tuple `(low, high, max_cons_low, max_cons_high)`: number of samples equal to `low` and `high`, and the length of the longest run of consecutive samples equal to `low` and to `high`.
"""
function saturation end
export saturation

function saturation(input::RadiationDetectorDSP.SamplesOrWaveform, low::RadiationDetectorDSP.RealQuantity, high::RadiationDetectorDSP.RealQuantity)
    X_axis, Y = RadiationDetectorDSP._get_axis_and_signal(input)
    from = firstindex(X_axis)
    until = lastindex(X_axis)
    _saturation_impl(X_axis, Y, from:until, low, high)
end

function saturation(input::RadiationDetectorDSP.SamplesOrWaveform, start::RadiationDetectorDSP.RealQuantity, stop::RadiationDetectorDSP.RealQuantity, low::RadiationDetectorDSP.RealQuantity, high::RadiationDetectorDSP.RealQuantity)
    X_axis, Y = RadiationDetectorDSP._get_axis_and_signal(input)
    # ToDo: Lower numerical precision of x-axis to y-axis, if x-axis is a range
    first_x, step_x = first(X_axis), step(X_axis)
    from = round(Int, ustrip(NoUnits, (start - first_x) / step_x)) + firstindex(X_axis)
    until = round(Int, ustrip(NoUnits, (stop - first_x) / step_x)) + firstindex(X_axis)
    _saturation_impl(X_axis, Y, from:until, low, high)
end

function _saturation_impl(X::AbstractArray{<:RadiationDetectorDSP.RealQuantity}, Y::AbstractArray{<:RadiationDetectorDSP.RealQuantity}, idxs::AbstractUnitRange{<:Integer}, low::RadiationDetectorDSP.RealQuantity, high::RadiationDetectorDSP.RealQuantity)
    @assert axes(X) == axes(Y)
    @assert firstindex(X) <= first(idxs) <= last(idxs) <= lastindex(X)
    @assert firstindex(Y) <= first(idxs) <= last(idxs) <= lastindex(Y)

    n_low       = 0
    n_high      = 0
    n_cons_low  = 0
    n_cons_high = 0
    n_const_low_counter = 0
    n_const_high_counter = 0
    @inbounds @fastmath @simd for i in idxs
        if Y[i] == low
            n_low += 1
            n_const_low_counter += 1
            n_cons_high = max(n_cons_high, n_const_high_counter)
            n_const_high_counter = 0
        elseif Y[i] == high
            n_high += 1
            n_const_high_counter += 1
            n_cons_low = max(n_cons_low, n_const_low_counter)
            n_const_low_counter = 0
        else
            n_cons_low = max(n_cons_low, n_const_low_counter)
            n_const_low_counter = 0
            n_cons_high = max(n_cons_high, n_const_high_counter)
            n_const_high_counter = 0
        end
    end
    n_cons_low = max(n_cons_low, n_const_low_counter)
    n_cons_high = max(n_cons_high, n_const_high_counter)
    (
        low = n_low,
        high = n_high,
        max_cons_low = n_cons_low,
        max_cons_high = n_cons_high
    )
end
