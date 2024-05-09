
"""
struct MultiIntersect <: Function
    Finds the x values at which a signal exceeds multiples of 10 percent of its 
    maximum value.
    Constructors:
    * ```$(FUNCTIONNAME)(; fields...)```
    Fields:
    $(TYPEDFIELDS)
"""
Base.@kwdef struct MultiIntersect{T<:RadiationDetectorDSP.RealQuantity} <: Function 
    "minimum time-over-threshold"
    mintot::T = 4
end

export MultiIntersect

function (f::MultiIntersect)(input::RadiationDetectorDSP.SamplesOrWaveform)

    X_axis, Y = RadiationDetectorDSP._get_axis_and_signal(input)
    tmp = f.mintot / step(X_axis)
    min_n_over_thresh = max(1, round(Int, ustrip(NoUnits, tmp)))
    _find_intersects_for_ann(X_axis, Y, min_n_over_thresh)
end

function _find_intersects_for_ann(X::AbstractVector{T}, 
    Y::AbstractVector{U}, min_n_over_thresh::Int
    ) where {T <: RealQuantity, U <: RealQuantity}
    
    thresholds = (1:9)/10 .* maximum(Y)
    l = length(thresholds)
    @assert axes(X) == axes(Y)

    # ToDo: What if eltype(Y) is based on ForwardDiff.Dual, but eltype(X) is not?
    R = float(T)
    intersect_x::Array{R} = zeros(R, l)


    isempty(Y) && return intersect_x

    # GPU-friendly branch-free code:

    cand_pos::Int = firstindex(Y) + 1
    intersect_pos::Array{Int} = fill(firstindex(Y), l)
    y_high_counter::Int = ifelse(first(Y) > first(thresholds), min_n_over_thresh + 1, 0)
    intersect_counter::Int = 1
    i::Int = firstindex(Y)
    @inbounds while (i <= lastindex(Y) && intersect_counter < l + 1)
        y = Y[i]
        y_is_high = y >= thresholds[intersect_counter]
        first_high_y = y_high_counter == 0
        cand_pos = ifelse(y_is_high && first_high_y, i, cand_pos)
        y_high_counter = ifelse(y_is_high, y_high_counter + 1, 0)
        new_intersect_found = y_high_counter == min_n_over_thresh
        pos = ifelse(new_intersect_found, cand_pos, intersect_pos[intersect_counter])
        intersect_pos[intersect_counter] = pos

        # reset running index i to previous found intersect pos
        i = ifelse(new_intersect_found, pos, i + 1)
        intersect_counter = ifelse(new_intersect_found, intersect_counter + 1, intersect_counter)
        y_high_counter = ifelse(new_intersect_found, 0, y_high_counter)
    end
    
    at_left_boundary = all(intersect_pos .== firstindex(Y))
    errmsg = "cannot interpolate intersect on left boundary"
    @assert !at_left_boundary errmsg

    for j in eachindex(intersect_pos)
        idx = intersect_pos[j]
        x_l = X[idx - 1]
        x_r = X[idx]
        y_l = Y[idx - 1]
        y_r = Y[idx]
        intersect_cand_x = R(thresholds[j] - y_l) * R(x_r - x_l) / R(y_r - y_l) + R(x_l)
        @assert y_l <= thresholds[j] <= y_r
        intersect_x[j] = intersect_cand_x
    end
    return intersect_x
end