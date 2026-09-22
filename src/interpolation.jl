# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    extrema3points(y1::T,y2::T,y3::T)::T where {T<:AbstractFloat}

Calculate the extrema of a parabola defined by three points.
"""
function extrema3points(y1::T,y2::T,y3::T)::T where {T<:RadiationDetectorDSP.RealQuantity}
    y1 - (y3-4y2+3y1)^2/(8*(y3-2y2+y1))
end

"""
    get_wvf_maximum(signal::AbstractSamples, start::Real, stop::Real)
    get_wvf_maximum(signal::RDWaveform, start::RealQuantity, stop::RealQuantity)

Get the maximum of a `signal` in the interval (`start`,`stop`) by quadaratic interpolation.
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
        if 1 < ind_max < length(idxs)
            wf_max = extrema3points(view(Y[idxs], ind_max-1:ind_max+1)...)
        else
            wf_max = Y[idxs][ind_max]
        end
    end

    return wf_max
end

"""
    get_wvf_maximum_time(signal::AbstractSamples, start::Real, stop::Real)
    get_wvf_maximum_time(signal::RDWaveform, start::RealQuantity, stop::RealQuantity)

Estimate the time of the maximum of `signal` in the interval (`start`, `stop`) by quadratic interpolation of the maximum sample and its two neighbors.
"""
function get_wvf_maximum_time end
export get_wvf_maximum_time

function get_wvf_maximum_time(input::RadiationDetectorDSP.SamplesOrWaveform, start::RadiationDetectorDSP.RealQuantity, stop::RadiationDetectorDSP.RealQuantity)
    X_axis, Y = RadiationDetectorDSP._get_axis_and_signal(input)
    first_x, step_x = first(X_axis), step(X_axis)
    from = round(Int, ustrip(NoUnits, (start - first_x) / step_x)) + firstindex(X_axis)
    until = round(Int, ustrip(NoUnits, (stop - first_x) / step_x)) + firstindex(X_axis)
    _get_wvf_maximum_time_impl(X_axis, Y, from:until)
end

function _get_wvf_maximum_time_impl(X::AbstractArray{<:RadiationDetectorDSP.RealQuantity}, Y::AbstractArray{<:RadiationDetectorDSP.RealQuantity}, idxs::AbstractUnitRange{<:Integer})
    @assert axes(X) == axes(Y)
    @assert firstindex(X) <= first(idxs) <= last(idxs) <= lastindex(X)
    @assert firstindex(Y) <= first(idxs) <= last(idxs) <= lastindex(Y)

    @inbounds begin
        ind_max = argmax(Y[idxs])
        maxidx = first(idxs) + ind_max - 1
        if 1 < ind_max < length(idxs)
            y1, y2, y3 = Y[maxidx-1], Y[maxidx], Y[maxidx+1]
            X[maxidx] + step(X) * (y1 - y3) / (2 * (y1 - 2y2 + y3))
        else
            X[maxidx]
        end
    end
end
