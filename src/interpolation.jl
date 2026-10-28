# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    extrema3points(y1::T, y2::T, y3::T)::T where {T<:RealQuantity}

Value at the extremum of the parabola through three equidistant samples `y1`, `y2` and `y3`.
"""
function extrema3points(y1::T,y2::T,y3::T)::T where {T<:RadiationDetectorDSP.RealQuantity}
    y1 - (y3-4y2+3y1)^2/(8*(y3-2y2+y1))
end

"""
    get_wvf_maximum(signal::SamplesOrWaveform, start::RealQuantity, stop::RealQuantity)

Maximum of `signal` between `start` and `stop` and its time, from a parabola through the maximum sample and its two neighbors.

The interval includes the samples nearest to `start` and `stop`. If the maximum sample is the first or last sample of the interval, its value and time are returned without interpolation.

# Returns
- named tuple `(max, t)`
"""
function get_wvf_maximum end
export get_wvf_maximum

function get_wvf_maximum(input::RadiationDetectorDSP.SamplesOrWaveform, start::RadiationDetectorDSP.RealQuantity, stop::RadiationDetectorDSP.RealQuantity)
    X_axis, Y = RadiationDetectorDSP._get_axis_and_signal(input)
    # ToDo: Lower numerical precision of x-axis to y-axis, if x-axis is a range
    first_x, step_x = first(X_axis), step(X_axis)
    from = round(Int, ustrip(NoUnits, (start - first_x) / step_x)) + firstindex(X_axis)
    until = round(Int, ustrip(NoUnits, (stop - first_x) / step_x)) + firstindex(X_axis)
    _get_wvf_maximum_impl(X_axis, Y, from:until)
end

function _get_wvf_maximum_impl(X::AbstractArray{<:RadiationDetectorDSP.RealQuantity}, Y::AbstractArray{<:RadiationDetectorDSP.RealQuantity}, idxs::AbstractUnitRange{<:Integer})
    @assert axes(X) == axes(Y)
    @assert firstindex(X) <= first(idxs) <= last(idxs) <= lastindex(X)
    @assert firstindex(Y) <= first(idxs) <= last(idxs) <= lastindex(Y)
    

    @inbounds begin
        ind_max = argmax(Y[idxs])
        maxidx = first(idxs) + ind_max - 1
        if 1 < ind_max < length(idxs)
            y1, y2, y3 = Y[maxidx-1], Y[maxidx], Y[maxidx+1]
            wvf_max = extrema3points(y1, y2, y3)
            δ = (y1 - y3) / (2 * (y1 - 2y2 + y3))
        else
            wvf_max = Y[maxidx]
            δ = zero(Y[maxidx] / Y[maxidx])
        end
        t_max = X[maxidx] + step(X) * δ
    end

    (max = wvf_max, t = t_max)
end
