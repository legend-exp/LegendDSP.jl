# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    struct IntersectMaximum <: Function
Finds the up-crossings of a signal `Y` over a threshold and picks the maximum in a
time window after each crossing.

Calling an `IntersectMaximum` on `(input, threshold)` returns a `NamedTuple`
`(x, x_tot, max, multiplicity)`: `x` the up-crossing time, `x_tot` the
time-over-threshold (down-crossing − up-crossing), `max` the maximum in the window,
and `multiplicity` the number of crossings found.

With `interpolation = true` (default) the up-/down-crossings are refined by linear
interpolation and the maximum by a 3-point parabola fit (`extrema3points`); with
`interpolation = false` the first sample over/under threshold and the highest raw
sample are used.

Constructors:
* ```$(FUNCTIONNAME)(; fields...)```
* ```$(FUNCTIONNAME)(mintot, maxtot; interpolation = true)```
Fields:
$(TYPEDFIELDS)
"""
Base.@kwdef struct IntersectMaximum{T<:RadiationDetectorDSP.RealQuantity} <: Function
    "minimum time-over-threshold"
    mintot::T = 4
    "maximum time-over-threshold for max to appear"
    maxtot::T = 100
    "if `true` (default), refine crossings by linear interpolation and the maximum by a 3-point parabola fit; if `false`, use the first sample over/under threshold and the highest raw sample"
    interpolation::Bool = true
end
export IntersectMaximum

# Positional `mintot`/`maxtot` with optional `interpolation` keyword: keeps existing
# `IntersectMaximum(mintot, maxtot)` calls working and lets `interpolation` be set
# without switching to all-keyword construction.
IntersectMaximum(mintot::T, maxtot::T; interpolation::Bool = true) where {T<:RadiationDetectorDSP.RealQuantity} =
    IntersectMaximum{T}(mintot, maxtot, interpolation)

function (f::IntersectMaximum)(input::RadiationDetectorDSP.SamplesOrWaveform, threshold::RadiationDetectorDSP.RealQuantity)
    X_axis, Y = RadiationDetectorDSP._get_axis_and_signal(input)
    min_n_over_thresh = max(1, round(Int, ustrip(NoUnits, f.mintot / step(X_axis))))
    max_n_for_max = max(1, round(Int, ustrip(NoUnits, f.maxtot / step(X_axis))))
    _find_intersect_maximum_impl(X_axis, Y, threshold, min_n_over_thresh, max_n_for_max, f.interpolation)
end
function _find_intersect_maximum_impl(X::AbstractVector{<:RadiationDetectorDSP.RealQuantity}, Y::AbstractVector{<:RadiationDetectorDSP.RealQuantity}, threshold::RadiationDetectorDSP.RealQuantity, min_n_over_thresh::Int, max_n_for_max::Int, interpolation::Bool)
    @assert axes(X) == axes(Y)
    # ToDo: What if eltype(Y) is based on ForwardDiff.Dual, but eltype(X) is not?
    R = float(eltype(X))
    M = float(eltype(Y))

    if isempty(Y)
        return (
            x = R[],
            x_tot = R[],
            max = M[],
            multiplicity = 0
        )
    end

    # Find up-crossings: positions where signal goes above threshold for at least min_n_over_thresh consecutive samples
    cand_pos::Int = firstindex(Y) + 1
    up_positions::Array{Int} = Int[]
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

    # Preallocate output arrays. The down-crossing is computed per pulse as a local
    # (`x_high_i`): it is needed for x_tot but is no longer part of the returned tuple.
    intersect_x = zeros(R, n)
    x_tot = zeros(R, n)
    intersect_max = zeros(M, n)

    for (i, up_pos) in enumerate(up_positions)
        # Up-crossing: linear interpolation (interpolation=true) or first sample ≥ threshold
        @inbounds begin
            if interpolation
                x_l = X[up_pos - 1]; x_r = X[up_pos]
                y_l = Y[up_pos - 1]; y_r = Y[up_pos]
                intersect_x[i] = R(threshold - y_l) * R(x_r - x_l) / R(y_r - y_l) + R(x_l)
            else
                intersect_x[i] = R(X[up_pos])
            end
        end

        # Find maximum in window after up-crossing
        from = max(up_pos - 2, firstindex(Y))
        until = min(up_pos + max_n_for_max, lastindex(Y))
        idxs = from:until
        @inbounds begin
            ind_max = argmax(Y[idxs])
            # 3-point parabola fit only when interpolating and the max is not at a window edge
            if interpolation && 1 < ind_max < length(idxs)
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

        # Down-crossing time (local; only used for x_tot): interpolated or first sample
        # < threshold; falls back to the waveform end if the signal never crosses back.
        if down_pos > firstindex(Y)
            @inbounds begin
                if interpolation
                    x_l = X[down_pos - 1]; x_r = X[down_pos]
                    y_l = Y[down_pos - 1]; y_r = Y[down_pos]
                    x_high_i = R(threshold - y_l) * R(x_r - x_l) / R(y_r - y_l) + R(x_l)
                else
                    x_high_i = R(X[down_pos])
                end
            end
        else
            x_high_i = R(last(X))
        end

        x_tot[i] = x_high_i - intersect_x[i]
    end

    return (
        x = intersect_x,
        x_tot = x_tot,
        max = intersect_max,
        multiplicity = n
    )
end