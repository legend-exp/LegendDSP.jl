# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    struct IntersectMaximum <: Function
Finds the intersects of a Y with a threshold and picking the maximum in a given time window.
Constructors:
* ```$(FUNCTIONNAME)(; fields...)```
Fields:
$(TYPEDFIELDS)
"""
Base.@kwdef struct IntersectMaximum{T<:RadiationDetectorDSP.RealQuantity} <: Function
    "minimum time-over-threshold"
    mintot::T = 4
    "maximum time-over-threshold for max to appear"
    maxtot::T = 100
end
export IntersectMaximum
function (f::IntersectMaximum)(input::RadiationDetectorDSP.SamplesOrWaveform, threshold::RadiationDetectorDSP.RealQuantity)
    X_axis, Y = RadiationDetectorDSP._get_axis_and_signal(input)
    min_n_over_thresh = max(1, round(Int, ustrip(NoUnits, f.mintot / step(X_axis))))
    max_n_for_max = max(1, round(Int, ustrip(NoUnits, f.maxtot / step(X_axis))))
    _find_intersect_maximum_impl(X_axis, Y, threshold, min_n_over_thresh, max_n_for_max)
end
function _find_intersect_maximum_impl(X::AbstractVector{<:RadiationDetectorDSP.RealQuantity}, Y::AbstractVector{<:RadiationDetectorDSP.RealQuantity}, threshold::RadiationDetectorDSP.RealQuantity, min_n_over_thresh::Int, max_n_for_max::Int)
    @assert axes(X) == axes(Y)
    # ToDo: What if eltype(Y) is based on ForwardDiff.Dual, but eltype(X) is not?
    R = float(eltype(X))
    if isempty(Y)
        return (
            x = R(NaN),
            multiplicity = -1
        )
    end
    # GPU-friendly branch-free code:
    cand_pos::Int = firstindex(Y) + 1
    intersect_pos_arr::Array{Int} = Int[]
    y_high_counter::Int = ifelse(first(Y) > threshold, min_n_over_thresh + 1, 0)
    n_intersects::Int = 0
    @inbounds for i in eachindex(Y)
        y = Y[i]
        y_is_high = y >= threshold
        first_high_y = y_high_counter == 0
        cand_pos = ifelse(y_is_high && first_high_y, i, cand_pos)
        y_high_counter = ifelse(y_is_high, y_high_counter + 1, 0)
        new_intersect_found = y_high_counter == min_n_over_thresh
        n_intersects = ifelse(new_intersect_found, n_intersects + 1, n_intersects)
        # no_previous_intersecs = n_intersects == 1
        # first_intersect_found = new_intersect_found && no_previous_intersecs
        # intersect_pos = ifelse(first_intersect_found, cand_pos, intersect_pos)
        if new_intersect_found
            push!(intersect_pos_arr, cand_pos)
        end
    end
    # @assert intersect_pos > firstindex(Y)
    intersect_cand_x = zeros(R, length(intersect_pos_arr))
    # Linear interpolation:
    for (i, intersect_pos) in enumerate(intersect_pos_arr)
        x_l = X[intersect_pos - 1]
        x_r = X[intersect_pos]
        y_l = Y[intersect_pos - 1]
        y_r = Y[intersect_pos]
        intersect_cand_x[i] = R(threshold - y_l) * R(x_r - x_l) / R(y_r - y_l) + R(x_l)
    end
    # @assert  y_l <= threshold <= y_r || n_intersects == 0
    # TODO: return NaN if no intersect found but make sure it is compatible with the other routines
    # intersect_x  = ifelse(n_intersects > 0, intersect_cand_x, zero(intersect_cand_x))
    n_intersects = ifelse(n_intersects > 0, n_intersects, 0)
    intersect_max = zeros(length(intersect_pos_arr))
    if n_intersects > 0
        for (i, intersect_pos) in enumerate(intersect_pos_arr)
            from = intersect_pos - 2
            until = intersect_pos + max_n_for_max
            idxs = from:until
            @inbounds begin
                ind_max = argmax(Y[idxs])
                if 2 < ind_max < length(X) - 1
                    intersect_max[i] = extrema3points(view(Y[idxs], ind_max-1:ind_max+1)...)
                else
                    intersect_max[i] = Y[ind_max]
                end
            end
        end
    end
    return (
        x = intersect_cand_x,
        max = intersect_max,
        multiplicity = n_intersects
    )
end