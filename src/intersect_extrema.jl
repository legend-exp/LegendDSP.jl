# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    struct IntersectExtrema{T<:RealQuantity} <: Function

Find the regions of a signal beyond a threshold, above `threshold` or below `-threshold`, and the extrema in each of them.

A region starts at the first sample beyond the threshold and ends at the first sample back inside it. It is counted if it lasts at least `mintot` and touches neither the first nor the last sample of the signal. Every region gives at least one extremum. Two extrema of a region are both kept if the signal between them comes back by more than `valley` from the smaller one; otherwise only the larger one is kept. The signal is expected to have its baseline at zero.

Call an instance as

    (f::IntersectExtrema)(signal::SamplesOrWaveform, threshold::RealQuantity, valley::RealQuantity)

It returns a named tuple of vectors with one element per extremum, sorted by time:
- `t`: time of the extremum, from a parabola through the extremal sample and its two neighbors
- `amp`: value of the signal at the extremum, from the same parabola, positive for a maximum and negative for a minimum
- `t_start`: time of the first sample of the region beyond the threshold
- `t_end`: time of the first sample after the region

Extrema of the same region have the same `t_start` and `t_end`. An `ArgumentError` is thrown if `threshold` or `valley` is negative or not a number.

`mintot` is converted to a number of samples using the step of the time axis. For a time axis with a unit it needs to be a quantity of the same dimension; the default of 4 is a plain number.

Constructors:

* ```$(FUNCTIONNAME)(; fields...)```
* ```$(FUNCTIONNAME)(mintot)```

Fields:

$(TYPEDFIELDS)
"""
Base.@kwdef struct IntersectExtrema{T<:RadiationDetectorDSP.RealQuantity} <: Function
    "minimum time beyond the threshold"
    mintot::T = 4
end
export IntersectExtrema

function (f::IntersectExtrema)(input::RadiationDetectorDSP.SamplesOrWaveform, threshold::RadiationDetectorDSP.RealQuantity, valley::RadiationDetectorDSP.RealQuantity)
    X_axis, Y = RadiationDetectorDSP._get_axis_and_signal(input)
    min_n_over_thresh = max(1, round(Int, ustrip(NoUnits, f.mintot / step(X_axis))))
    _find_intersect_extrema_impl(X_axis, Y, threshold, valley, min_n_over_thresh)
end

function _find_intersect_extrema_impl(X::AbstractVector{<:RadiationDetectorDSP.RealQuantity}, Y::AbstractVector{<:RadiationDetectorDSP.RealQuantity}, threshold::RadiationDetectorDSP.RealQuantity, valley::RadiationDetectorDSP.RealQuantity, min_n_over_thresh::Int)
    @assert axes(X) == axes(Y)
    (threshold >= zero(threshold) && valley >= zero(valley)) || throw(ArgumentError("IntersectExtrema needs a non-negative threshold and valley, got $threshold and $valley"))
    R = float(eltype(X))
    M = float(eltype(Y))
    t = R[]
    amp = M[]
    t_start = R[]
    t_end = R[]

    @inbounds for s in (1, -1)
        i = firstindex(Y) + 1
        while i <= lastindex(Y)
            if !(s * Y[i] > threshold && s * Y[i-1] <= threshold)
                i += 1
                continue
            end
            from = i
            while i <= lastindex(Y) && s * Y[i] > threshold
                i += 1
            end
            (i <= lastindex(Y) && i - from >= min_n_over_thresh) || continue

            extrema_pos = Int[]
            for c in from:i-1
                (s * Y[c] > s * Y[c-1] && s * Y[c] >= s * Y[c+1]) || continue
                if isempty(extrema_pos) || min(s * Y[last(extrema_pos)], s * Y[c]) - minimum(k -> s * Y[k], last(extrema_pos):c) > valley
                    push!(extrema_pos, c)
                elseif s * Y[c] > s * Y[last(extrema_pos)]
                    extrema_pos[end] = c
                end
            end

            for c in extrema_pos
                y1, y2, y3 = Y[c-1], Y[c], Y[c+1]
                push!(t, X[c] + step(X) * (y1 - y3) / (2 * (y1 - 2y2 + y3)))
                push!(amp, extrema3points(y1, y2, y3))
                push!(t_start, X[from])
                push!(t_end, X[i])
            end
        end
    end

    o = sortperm(t)
    return (
        t = t[o],
        amp = amp[o],
        t_start = t_start[o],
        t_end = t_end[o]
    )
end
