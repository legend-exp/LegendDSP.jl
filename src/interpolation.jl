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
        if 2 < ind_max < length(X) - 1
            wf_max = extrema3points(view(Y, ind_max-1:ind_max+1)...)
        else
            wf_max = Y[ind_max]
        end
    end

    return wf_max
end
