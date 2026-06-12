# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    histogram_baseline(signal::AbstractSamples; kwargs...)
    histogram_baseline(signal::RDWaveform; kwargs...)

Estimate the baseline and noise of a `signal` from the y-projection histogram
of its samples. The mode (most common value) is the baseline; the half width
at half maximum (HWHM) of the peak gives the noise.

# Binning (two modes)
- **Fixed bin width** — pass `bin_width` (in signal/ADC units). Bins are laid on
  a grid whose *centres are multiples of `bin_width`*, so `bin_width = 1` puts
  every integer ADC value on its own bin centre (…, 2999, 3000, 3001, …) — the
  natural choice for an integer-quantised DAQ. The grid spans the full data
  range; *all* samples are used.
- **Robust, fixed bin count** (default, `bin_width = nothing`) — the histogram
  range is `median ± n_sigma · σ_robust`, where `σ_robust = 1.4826 · MAD`
  (median absolute deviation). This auto-zooms onto the baseline peak and is
  *insensitive to extreme outliers* (a saturated sample neither moves the
  median nor the MAD), so the `nbins` bins stay narrow in ADC units. Note that
  a plain standard deviation would NOT work here — one outlier inflates it and
  the window blows up again, which is exactly why MAD is used.

The counts are (optionally) smoothed with a moving-average kernel to suppress
Poisson noise before the mode/HWHM are extracted. The asymmetry
`hwhm_right - hwhm_left` flags pulse contamination (pulses widen the right
flank).

# Keyword arguments
- `nbins::Integer = 100`: number of bins (robust mode only; ignored if
  `bin_width` is given).
- `bin_width = nothing`: if set (signal/ADC units), use fixed-width bins over
  the full data range instead of the robust `nbins` window.
- `n_sigma::Real = 5.0`: half-width of the robust histogram window in units of
  `σ_robust` (robust mode only).
- `smooth::Integer = 1`: width of the centered moving-average kernel (positive
  odd integer; `1` disables smoothing — the default, since an integer-ADC
  histogram needs none). Raise it for sparse/continuous data if the peak is noisy.
- `interpolation::Bool = true`: if `true`, the mode is parabola-refined to
  sub-bin precision and the HWHM crossings are linearly interpolated. If
  `false`, the mode is just the highest (smoothed) bin and each HWHM is the
  distance to the first bin that drops below half-maximum — simpler and faster.
- `min`, `max`: hard pre-filter; only samples inside `[min, max]` are used at
  all. Defaults to `±Inf`, i.e. no filtering.
- `return_hist::Bool = false`: if `true`, add a `hist` field exposing the exact
  binned histogram (over the window actually used) for plotting. No extra
  allocation when `false`.

# Returns
`NamedTuple` with:
- `mode`         — baseline estimate (parabola-refined if `interpolation`)
- `hwhm_left`    — distance from mode to half-maximum crossing on the left
- `hwhm_right`   — distance from mode to half-maximum crossing on the right
- `sigma_left`   — `hwhm_left / √(2 ln 2)` (Gaussian-equivalent σ, left side)
- `sigma_right`  — `hwhm_right / √(2 ln 2)` (Gaussian-equivalent σ, right side)
- `sigma`        — average of `sigma_left` and `sigma_right`, a drop-in
                   replacement for `signalstats(...).sigma`
- `mode_count`   — raw (unsmoothed) sample count in the mode bin
- `n_used`       — number of samples that passed the `[min, max]` filter
- `hist`         — *only if `return_hist`*: a `NamedTuple`
                   `(centers, counts, smoothed, lo, hi, bin_width, mode_bin, half_max)`
                   describing the binned histogram. `centers` are bin centres,
                   `counts`/`smoothed` the raw/smoothed bin contents, `lo`/`hi`
                   the window binned, and `half_max = smoothed[mode_bin]/2`.
"""
function histogram_baseline end
export histogram_baseline

function histogram_baseline(input::RadiationDetectorDSP.SamplesOrWaveform{T};
        nbins::Integer = 100,
        bin_width::Union{Nothing,RadiationDetectorDSP.RealQuantity} = nothing,
        n_sigma::Real = 5.0,
        smooth::Integer = 1,
        interpolation::Bool = true,
        min::RadiationDetectorDSP.RealQuantity = -Inf * unit(T),
        max::RadiationDetectorDSP.RealQuantity = Inf * unit(T),
        return_hist::Bool = false,
    ) where T
    _, Y = RadiationDetectorDSP._get_axis_and_signal(input)
    bw = isnothing(bin_width) ? nothing : T(bin_width)
    _histogram_baseline_impl(Y, T(min), T(max), Int(nbins), Int(smooth), bw,
                             float(n_sigma), interpolation, return_hist)
end

function _histogram_baseline_impl(Y::AbstractArray{T}, _min::T, _max::T, nbins::Int,
        smooth::Int, bin_width::Union{Nothing,T}, n_sigma::Real, interpolation::Bool,
        return_hist::Bool) where {T <: RealQuantity}
    @assert nbins >= 2 "nbins must be at least 2"
    @assert smooth >= 1 && isodd(smooth) "smooth must be a positive odd integer"
    @assert n_sigma > 0 "n_sigma must be positive"

    # First pass: extent and count of in-range samples.
    y_lo = typemax(T)
    y_hi = typemin(T)
    n_used = 0
    @inbounds for i in eachindex(Y)
        y = Y[i]
        if _min <= y <= _max
            n_used += 1
            y_lo = ifelse(y < y_lo, y, y_lo)
            y_hi = ifelse(y > y_hi, y, y_hi)
        end
    end

    if n_used == 0 || !(y_hi > y_lo)
        bl = (n_used == 0 ? zero(T) : y_lo)
        base = (mode = bl, hwhm_left = zero(T), hwhm_right = zero(T),
                sigma_left = zero(T), sigma_right = zero(T), sigma = zero(T),
                mode_count = n_used, n_used = n_used)
        return_hist || return base
        # Degenerate input: emit an empty histogram so the `hist` field set is
        # stable regardless of input (downstream plots read `r.hist` unconditionally).
        return merge(base, (hist = (centers = T[], counts = Int[], smoothed = Float64[],
                                    lo = bl, hi = bl, bin_width = zero(T),
                                    mode_bin = 0, half_max = 0.0),))
    end

    # Decide histogram window [lo, hi] and bin width.
    if bin_width === nothing
        # Robust window: median ± n_sigma · (1.4826·MAD). Both the centre
        # (median) and the scale (MAD) ignore extreme outliers, so the window
        # hugs the baseline peak regardless of saturated samples / big pulses.
        c, s = _robust_center_scale(Y, _min, _max, n_used)
        lo = max(c - n_sigma * s, _min)
        hi = min(c + n_sigma * s, _max)
        if !(hi > lo)
            lo, hi = y_lo, y_hi   # MAD == 0 (flat) → fall back to data extent
        end
        nb = nbins
        bw = (hi - lo) / nb
    else
        # Fixed bin width, snapped to a grid whose centres are multiples of bw:
        # bin j is centred on (k_lo+j-1)·bw, so bw = 1 ADC puts every integer value
        # on its own centre (…, 3000, 3001, …). The DAQ digitises to whole ADC, so
        # each quantised level gets exactly one bin. All in-range samples are used.
        bw = bin_width
        k_lo = round(Int, ustrip(NoUnits, y_lo / bw))
        k_hi = round(Int, ustrip(NoUnits, y_hi / bw))
        lo = (k_lo - 0.5) * bw      # left edge = first centre − bw/2
        hi = (k_hi + 0.5) * bw      # right edge = last centre + bw/2
        nb = clamp(k_hi - k_lo + 1, 2, 10_000_000)
    end

    # Build histogram over [lo, hi]; samples outside the window are dropped
    # (not clamped into edge bins, which would pile pulses onto a flank).
    counts = zeros(Int, nb)
    @inbounds for i in eachindex(Y)
        y = Y[i]
        if _min <= y <= _max && lo <= y <= hi
            idx = floor(Int, ustrip(NoUnits, (y - lo) / bw)) + 1
            idx = clamp(idx, 1, nb)
            counts[idx] += 1
        end
    end

    smoothed = _moving_average(counts, smooth)
    mode_bin = argmax(smoothed)

    # Mode position. With interpolation, refine sub-bin via a parabola through
    # the three central bins; vertex offset δ = 0.5·(c₋ - c₊)/(c₋ - 2c₀ + c₊).
    δ_bin = 0.0
    if interpolation && 1 < mode_bin < nb
        c_lo = smoothed[mode_bin - 1]
        c_0  = smoothed[mode_bin]
        c_hi = smoothed[mode_bin + 1]
        denom = c_lo - 2 * c_0 + c_hi
        if denom < 0
            δ_bin = 0.5 * (c_lo - c_hi) / denom
            δ_bin = clamp(δ_bin, -0.5, 0.5)
        end
    end
    mode_count = counts[mode_bin]
    mode = lo + (mode_bin - 0.5 + δ_bin) * bw

    half_max = smoothed[mode_bin] / 2
    if interpolation
        hwhm_left  = _half_max_crossing_left(smoothed,  mode_bin, half_max, bw, δ_bin)
        hwhm_right = _half_max_crossing_right(smoothed, mode_bin, half_max, bw, δ_bin)
    else
        hwhm_left  = _nearest_below_left(smoothed,  mode_bin, half_max, bw)
        hwhm_right = _nearest_below_right(smoothed, mode_bin, half_max, bw)
    end

    inv_kfwhm = oftype(float(one(T)) / one(T), 1 / sqrt(2 * log(2)))
    σ_l = hwhm_left  * inv_kfwhm
    σ_r = hwhm_right * inv_kfwhm
    σ   = (σ_l + σ_r) / 2

    base = (mode = mode, hwhm_left = hwhm_left, hwhm_right = hwhm_right,
            sigma_left = σ_l, sigma_right = σ_r, sigma = σ,
            mode_count = mode_count, n_used = n_used)
    return_hist || return base   # fast path: no extra allocation when not requested

    # Expose the exact binned histogram (over the window actually used) so a plot
    # can reproduce what the algorithm saw. Bin CENTRES (not edges) keep the mode
    # and HWHM markers on the same axis as the bars.
    centers = [lo + (i - 0.5) * bw for i in 1:nb]
    merge(base, (hist = (centers = centers, counts = counts, smoothed = smoothed,
                         lo = lo, hi = hi, bin_width = bw,
                         mode_bin = mode_bin, half_max = half_max),))
end

# Robust centre (median) and scale (1.4826·MAD) of the in-range samples.
# Copies the in-range values into a scratch buffer so the two `median!` calls
# (which reorder their input) don't touch the caller's data.
function _robust_center_scale(Y::AbstractArray{T}, _min::T, _max::T, n_used::Int) where {T <: RealQuantity}
    buf = Vector{T}(undef, n_used)
    k = 0
    @inbounds for i in eachindex(Y)
        y = Y[i]
        if _min <= y <= _max
            k += 1
            buf[k] = y
        end
    end
    c = median!(buf)
    @inbounds for i in eachindex(buf)
        buf[i] = abs(buf[i] - c)
    end
    kmad = oftype(float(one(T)), 1.4826)
    s = kmad * median!(buf)
    return c, s
end

# --- HWHM, interpolated (sub-bin) ---------------------------------------------
function _half_max_crossing_left(counts::AbstractVector{<:Real}, mode_bin::Int, half_max::Real, bin_width::T, δ_bin::Real) where {T}
    @inbounds for i in (mode_bin-1):-1:1
        if counts[i] < half_max
            c_lo = counts[i]; c_hi = counts[i+1]
            denom = c_hi - c_lo
            f = denom > 0 ? (half_max - c_lo) / denom : 0.5
            return ((mode_bin - 0.5 + δ_bin) - (i - 0.5 + f)) * bin_width
        end
    end
    return (mode_bin - 1 + δ_bin) * bin_width
end

function _half_max_crossing_right(counts::AbstractVector{<:Real}, mode_bin::Int, half_max::Real, bin_width::T, δ_bin::Real) where {T}
    n = length(counts)
    @inbounds for i in (mode_bin+1):n
        if counts[i] < half_max
            c_hi = counts[i-1]; c_lo = counts[i]
            denom = c_hi - c_lo
            f = denom > 0 ? (c_hi - half_max) / denom : 0.5
            return ((i - 1.5 + f) - (mode_bin - 0.5 + δ_bin)) * bin_width
        end
    end
    return (n - mode_bin + 0.5 - δ_bin) * bin_width
end

# --- HWHM, nearest-bin (no interpolation) -------------------------------------
# Distance from the mode bin to the first bin that drops below half-maximum.
function _nearest_below_left(counts::AbstractVector{<:Real}, mode_bin::Int, half_max::Real, bin_width::T) where {T}
    @inbounds for i in (mode_bin-1):-1:1
        counts[i] < half_max && return (mode_bin - i) * bin_width
    end
    return (mode_bin - 1) * bin_width
end

function _nearest_below_right(counts::AbstractVector{<:Real}, mode_bin::Int, half_max::Real, bin_width::T) where {T}
    n = length(counts)
    @inbounds for i in (mode_bin+1):n
        counts[i] < half_max && return (i - mode_bin) * bin_width
    end
    return (n - mode_bin) * bin_width
end

# Centered moving-average smoothing, returning Float64. Width must be odd ≥ 1.
function _moving_average(counts::AbstractVector{Int}, width::Int)
    n = length(counts)
    if width == 1
        return Float64.(counts)
    end
    half = (width - 1) ÷ 2
    out = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        s = 0
        cnt = 0
        for j in (i - half):(i + half)
            if 1 <= j <= n
                s += counts[j]
                cnt += 1
            end
        end
        out[i] = s / cnt
    end
    return out
end
