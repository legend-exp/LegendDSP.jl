# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    iterative_baseline(signal::AbstractSamples; kwargs...)
    iterative_baseline(signal::RDWaveform; kwargs...)

Iteratively compute the mean and σ of a `signal`, in each iteration restricting
the accepted samples to a ±`n_sigma`·σ band around the previous mean.

This is iterative σ-clipping (a standard robust-statistics technique): start
with a plain mean/σ over all samples, then iterate
```
sigma_new = stddev of samples in [mean_old - n_sigma·sigma_old, mean_old + n_sigma·sigma_old]
mean_new  = mean    of samples in the same window
```
until either `|mean_new - mean_old| < sigma_new / √(n_accepted - 1)` (the
mean's own standard error) or `max_iter` is reached. The first pass uses
all samples (no window restriction yet).

Useful when the baseline window contains a few pulses that would bias a
plain `signalstats` mean — the iterative cut converges to the underlying
noise distribution and rejects the pulse samples.

# Keyword arguments
- `n_sigma::Real = 3.0`  : half-width of the acceptance window in σ.
- `max_iter::Integer = 12`: maximum number of iterations.
- `return_mask::Bool = false`: if `true`, add an `accepted::BitVector` flagging
  the samples inside the final acceptance window that defined the returned
  mean/sigma. No extra allocation when `false`.

# Returns
`NamedTuple` with:
- `mean`        — final baseline estimate
- `sigma`       — final noise estimate (std of accepted samples)
- `n_iter`      — number of completed iterations
- `n_accepted`  — number of samples inside the final window
- `converged`   — `true` if the convergence criterion was met before `max_iter`
- `accepted`    — *only if `return_mask`*: a `BitVector` (length = #samples)
                  marking the samples inside the defining window;
                  `count(accepted) == n_accepted`.
"""
function iterative_baseline end
export iterative_baseline

function iterative_baseline(input::RadiationDetectorDSP.SamplesOrWaveform{T};
        n_sigma::Real = 3.0,
        max_iter::Integer = 12,
        return_mask::Bool = false,
    ) where T
    _, Y = RadiationDetectorDSP._get_axis_and_signal(input)
    _iterative_baseline_impl(Y, float(n_sigma), Int(max_iter), return_mask)
end

function _iterative_baseline_impl(Y::AbstractArray{T}, n_sigma::Real, max_iter::Int,
        return_mask::Bool) where {T <: RealQuantity}
    @assert max_iter >= 1 "max_iter must be at least 1"
    @assert n_sigma > 0  "n_sigma must be positive"

    N = length(Y)
    if N == 0
        base = (
            mean       = zero(T),
            sigma      = zero(T),
            n_iter     = 0,
            n_accepted = 0,
            converged  = false,
        )
        return return_mask ? merge(base, (accepted = falses(0),)) : base
    end

    # Initial pass: plain mean / sigma over all samples
    mean_old, sigma_old, n_acc = _windowed_meanstd(Y, typemin(T), typemax(T))
    if n_acc <= 1 || iszero(sigma_old)
        # Cannot iterate without a non-zero spread; return the plain estimate.
        base = (
            mean       = mean_old,
            sigma      = sigma_old,
            n_iter     = 1,
            n_accepted = n_acc,
            converged  = true,
        )
        # The returned estimate used ALL samples (no window) → everything accepted.
        return return_mask ? merge(base, (accepted = trues(N),)) : base
    end

    converged = false
    n_iter = 1
    mean_new = mean_old
    sigma_new = sigma_old
    # Track the window that defines the returned mean/sigma: it is the [lo, hi]
    # computed at the TOP of the final iteration, captured before mean_old/sigma_old
    # are overwritten at the loop's end.
    lo_final = mean_old - n_sigma * sigma_old
    hi_final = mean_old + n_sigma * sigma_old
    while n_iter < max_iter
        lo = mean_old - n_sigma * sigma_old
        hi = mean_old + n_sigma * sigma_old
        lo_final, hi_final = lo, hi
        mean_new, sigma_new, n_acc = _windowed_meanstd(Y, lo, hi)
        n_iter += 1

        if n_acc <= 1 || iszero(sigma_new)
            # Window collapsed (all samples rejected or all identical) — stop.
            converged = true
            break
        end

        # Convergence: change in mean below the mean's own standard error
        tol = sigma_new / sqrt(oftype(float(one(sigma_new) / one(sigma_new)), n_acc - 1))
        if abs(mean_new - mean_old) < tol
            converged = true
            break
        end

        mean_old = mean_new
        sigma_old = sigma_new
    end

    base = (
        mean       = mean_new,
        sigma      = sigma_new,
        n_iter     = n_iter,
        n_accepted = n_acc,
        converged  = converged,
    )
    return_mask || return base   # fast path: no extra allocation when not requested

    # Reconstruct the accepted mask from the defining window. Same predicate and
    # bounds as `_windowed_meanstd`, so `count(accepted) == n_accepted` by construction.
    accepted = BitVector(undef, N)
    @inbounds @simd for i in eachindex(Y)
        accepted[i] = lo_final <= Y[i] <= hi_final
    end
    merge(base, (accepted = accepted,))
end

# Compute (mean, std, n) of samples in `Y` that fall inside `[lo, hi]`.
# Returns std as a population standard deviation (divides by `n_acc`).
function _windowed_meanstd(Y::AbstractArray{T}, lo::T, hi::T) where {T <: RealQuantity}
    sum_Y::float(typeof(zero(T))) = zero(float(T))
    sum_Y_sqr::float(typeof(zero(T) * zero(T))) = zero(float(T)) * zero(float(T))
    n_acc = 0
    @inbounds @fastmath @simd for i in eachindex(Y)
        y = Y[i]
        _include = lo <= y <= hi
        # Branchless accumulation: zero-out contributions when outside the window.
        y_in = ifelse(_include, y, zero(T))
        sum_Y = sum_Y + y_in
        sum_Y_sqr = fma(y_in, y_in, sum_Y_sqr)
        n_acc += _include
    end
    if n_acc == 0
        return zero(T), zero(T), 0
    end
    inv_n = inv(n_acc)
    mean_Y = sum_Y * inv_n
    var_Y = sum_Y_sqr * inv_n - mean_Y * mean_Y
    var_Y = ifelse(var_Y < zero(var_Y), zero(var_Y), var_Y)
    return mean_Y, sqrt(var_Y), n_acc
end
