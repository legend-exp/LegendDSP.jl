# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    extrema_over_threshold(input::SamplesOrWaveform, threshold::RealQuantity)

Find the sites of a SiPM trace: every local maximum above `threshold` and every local minimum below `-threshold`.
A site gets its time (parabolic refinement of the extremum), its signed extremum value and the time over threshold
of the excursion it lies in (both crossings linearly interpolated). A flat top counts once, at its first sample, and
only if the trace falls after it. Sites are sorted by time. `threshold` must not be negative.

Returns `(x, x_tot, max)`.
"""
function extrema_over_threshold(input::RadiationDetectorDSP.SamplesOrWaveform, threshold::RadiationDetectorDSP.RealQuantity)
    threshold >= zero(threshold) || throw(ArgumentError("extrema_over_threshold needs a non-negative threshold, got $threshold"))
    X, Y = RadiationDetectorDSP._get_axis_and_signal(input)
    R = float(eltype(X))
    M = float(eltype(Y))
    x = R[]
    x_tot = R[]
    site_max = M[]
    n = length(Y)
    @inbounds for s in (1, -1), i in 2:n-1
        v = s * Y[i]
        (v > threshold && v > s * Y[i-1] && v >= s * Y[i+1]) || continue
        # a flat top: the trace must fall after it, otherwise it is a step on an edge
        j = i
        while j < n && Y[j+1] == Y[i]
            j += 1
        end
        (j == n || s * Y[j+1] < v) || continue
        # the excursion above threshold around the extremum and its interpolated crossings
        a = i
        while a > 1 && s * Y[a-1] > threshold
            a -= 1
        end
        b = i
        while b < n && s * Y[b+1] > threshold
            b += 1
        end
        thr = s * threshold
        x_lo = a > 1 ? R(X[a-1] + (thr - Y[a-1]) * (X[a] - X[a-1]) / (Y[a] - Y[a-1])) : R(X[a])
        x_hi = b < n ? R(X[b] + (thr - Y[b]) * (X[b+1] - X[b]) / (Y[b+1] - Y[b])) : R(X[b])
        # parabola through the extremum and its neighbours
        y1, y2, y3 = Y[i-1], Y[i], Y[i+1]
        d = y1 - 2y2 + y3
        δ = d == 0 ? zero(M) : (y1 - y3) / (2d)
        push!(x, R(X[i] + δ * (X[i+1] - X[i])))
        push!(x_tot, x_hi - x_lo)
        push!(site_max, d == 0 ? M(y2) : M(y2 - (y3 - y1)^2 / (8d)))
    end
    o = sortperm(x)
    (x = x[o], x_tot = x_tot[o], max = site_max[o])
end
export extrema_over_threshold

"""
    struct SiteSolver <: Function

Solve the photon amplitudes at the sites of a SiPM trace on a long shaping from the response of one photon: the trace
read at the sites is the sum of the responses of all sites, `b = K A` with `K[i, j] = response(x_i - x_j)`, so
`A = K \\ b` removes what the neighbours contribute at every site. Same-sign sites closer than `merge_distance`
are merged first (centroid weighted with the read values; the long shaping cannot separate them), the solve is
sign constrained (the site whose amplitude contradicts its sign the most is set to zero and the rest solved
again) and a singular system falls back to the read values.

Called as `solver(y, pos, tot, sgn)` with the long trace `y`, the site positions `pos` (samples, 1-based), their
times over threshold `tot` and signs `sgn`; returns `(pos, tot, sgn, A, readoff)` after merging, sorted by position.
Every position must lie within the trace.

Fields:

$(TYPEDFIELDS)
"""
Base.@kwdef struct SiteSolver{TV <: AbstractVector{<:Real}} <: Function
    "response of one photon on the long shaping, sampled on the trace grid, peak 1"
    response::TV
    "index of the response peak"
    i_peak::Int
    "same-sign sites closer than this (samples) merge"
    merge_distance::Float64
end
export SiteSolver

function (f::SiteSolver)(y::AbstractVector{<:Real}, pos::AbstractVector{<:Real}, tot::AbstractVector{<:Real}, sgn::AbstractVector{<:Integer})
    n = length(y)
    all(p -> 1 <= round(Int, p) <= n, pos) || throw(ArgumentError("SiteSolver: site positions $(extrema(pos)) outside the trace of $n samples"))
    pos = Float64.(pos)
    tot = Float64.(tot)
    sgn = Int.(sgn)
    read(p) = Float64(y[round(Int, p)])
    b = [read(p) for p in pos]
    # merge the closest same-sign pair within merge_distance until none is left
    while length(pos) > 1
        ic = jc = 0
        dmin = f.merge_distance
        for k in 2:length(pos), i in 1:k-1
            sgn[i] == sgn[k] || continue
            d = abs(pos[i] - pos[k])
            d < dmin && (dmin = d; ic = i; jc = k)
        end
        ic == 0 && break
        wi, wj = abs(b[ic]), abs(b[jc])
        pos[ic] = wi + wj > 0 ? (wi * pos[ic] + wj * pos[jc]) / (wi + wj) : (pos[ic] + pos[jc]) / 2
        tot[ic] = max(tot[ic], tot[jc])
        b[ic] = read(pos[ic])
        deleteat!(pos, jc)
        deleteat!(tot, jc)
        deleteat!(sgn, jc)
        deleteat!(b, jc)
    end
    # sign-constrained solve
    kern(d) = (k = f.i_peak + round(Int, d); 1 <= k <= length(f.response) ? Float64(f.response[k]) : 0.0)
    A = zeros(length(pos))
    act = collect(eachindex(pos))
    while !isempty(act)
        K = [kern(pos[i] - pos[k]) for i in act, k in act]
        F = lu(K; check = false)
        Aa = issuccess(F) ? F \ b[act] : b[act]
        viol = [q for q in eachindex(act) if Aa[q] * sgn[act[q]] < 0]
        if isempty(viol)
            A[act] .= Aa
            break
        end
        deleteat!(act, viol[argmax(abs.(Aa[viol]))])
    end
    o = sortperm(pos)
    (pos = pos[o], tot = tot[o], sgn = sgn[o], A = A[o], readoff = b[o])
end
