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
            x = R[],
            x_high = R[],
            x_tot = R[],
            max = Float64[],
            multiplicity = 0
        )
    end

    # Find up-crossings: positions where signal goes above threshold for at least min_n_over_thresh consecutive samples
    cand_pos::Int = firstindex(Y) + 1
    up_positions = Int[]
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
        if new_intersect_found && cand_pos > firstindex(Y)
            push!(up_positions, cand_pos)
        end
    end

    n_intersects = max(n_intersects, 0)
    n = length(up_positions)

    # Preallocate output arrays
    intersect_x = zeros(R, n)
    x_high = zeros(R, n)
    x_tot = zeros(R, n)
    intersect_max = zeros(n)

    for (i, up_pos) in enumerate(up_positions)
        # Linear interpolation for up-crossing
        @inbounds begin
            x_l = X[up_pos - 1]; x_r = X[up_pos]
            y_l = Y[up_pos - 1]; y_r = Y[up_pos]
            intersect_x[i] = R(threshold - y_l) * R(x_r - x_l) / R(y_r - y_l) + R(x_l)
        end

        # Find maximum in window after up-crossing
        from = max(up_pos - 2, firstindex(Y))
        until = min(up_pos + max_n_for_max, lastindex(Y))
        idxs = from:until
        @inbounds begin
            ind_max = argmax(Y[idxs])
            if 1 < ind_max < length(idxs)
                intersect_max[i] = extrema3points(view(Y[idxs], ind_max-1:ind_max+1)...)
            else
                intersect_max[i] = Y[idxs][ind_max]
            end
        end

        # Find down-crossing: first sample below threshold after the confirmed up-crossing region
        down_pos = 0
        @inbounds for j in (up_pos + min_n_over_thresh):lastindex(Y)
            if Y[j] < threshold
                down_pos = j
                break
            end
        end

        # Linear interpolation for down-crossing (x_high)
        if down_pos > firstindex(Y)
            @inbounds begin
                x_l = X[down_pos - 1]; x_r = X[down_pos]
                y_l = Y[down_pos - 1]; y_r = Y[down_pos]
                x_high[i] = R(threshold - y_l) * R(x_r - x_l) / R(y_r - y_l) + R(x_l)
            end
        else
            # No down-crossing found: signal stays above threshold until end of waveform
            x_high[i] = R(last(X))
        end

        x_tot[i] = x_high[i] - intersect_x[i]
    end

    return (
        x = intersect_x,
        x_high = x_high,
        x_tot = x_tot,
        max = intersect_max,
        multiplicity = n_intersects
    )
end