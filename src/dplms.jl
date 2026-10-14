# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

# DPLMS optimum filter synthesis: waveform alignment, noise and signal matrices, solver.

# Index of the sample of the range `t` closest to `x`.
_sample_index(t::AbstractRange, x) = round(Int, ustrip(NoUnits, (x - first(t)) / step(t))) + firstindex(t)

"""
    align_waveform(wvf::RDWaveform, t_align, t_target)

Shift `wvf` in time so that the feature at `t_align` lands at `t_target`. The output has the
time axis and length of the input; its value at time `t` is the input at
`t + (t_align - t_target)`, interpolated with a cubic Catmull-Rom spline through the four
nearest samples. Positions outside the record repeat the first or last sample.
"""
function align_waveform(wvf::RDWaveform, t_align, t_target)
    shift = ustrip(NoUnits, (t_align - t_target) / step(wvf.time))
    RDWaveform(wvf.time, _shift_samples(wvf.signal, shift))
end
export align_waveform

"""
    align_waveforms(wvfs::ArrayOfRDWaveforms, t_align::AbstractVector, t_target)

[`align_waveform`](@ref) for every waveform with its own alignment time.
"""
function align_waveforms(wvfs::ArrayOfRDWaveforms, t_align::AbstractVector, t_target)
    axes(wvfs) == axes(t_align) || throw(DimensionMismatch("waveforms $(axes(wvfs)) and alignment times $(axes(t_align)) differ"))
    ArrayOfRDWaveforms([align_waveform(wvfs[i], t_align[i], t_target) for i in eachindex(wvfs, t_align)])
end
export align_waveforms

# Samples of `y` read `shift` samples later: out[i] = y(i + shift) by Catmull-Rom
# interpolation, edge samples repeated outside the record.
function _shift_samples(y::AbstractVector, shift::Real)
    out = similar(y, float(eltype(y)))
    i_first, i_last = firstindex(y), lastindex(y)
    k = floor(Int, shift)
    f = shift - k
    w0 = (-f^3 + 2f^2 - f) / 2
    w1 = (3f^3 - 5f^2 + 2) / 2
    w2 = (-3f^3 + 4f^2 + f) / 2
    w3 = (f^3 - f^2) / 2
    for i in eachindex(y)
        j = i + k
        p0 = y[clamp(j - 1, i_first, i_last)]
        p1 = y[clamp(j, i_first, i_last)]
        p2 = y[clamp(j + 1, i_first, i_last)]
        p3 = y[clamp(j + 2, i_first, i_last)]
        out[i] = w0 * p0 + w1 * p1 + w2 * p2 + w3 * p3
    end
    out
end

# Time axis shared by all waveforms of the batch; throws if any waveform has a different one.
function _shared_time_axis(wvfs::ArrayOfRDWaveforms)
    t = first(wvfs.time)
    all(==(t), wvfs.time) || throw(ArgumentError("waveforms do not share one time axis"))
    t
end

"""
    aligned_windows(wvfs::ArrayOfRDWaveforms, t_target, n_pre::Integer, n_post::Integer)

Samples from `n_pre` before to `n_post` after the sample at `t_target`, one column per
waveform. Throws if the window leaves the record, or if the waveforms do not share one
time axis.
"""
function aligned_windows(wvfs::ArrayOfRDWaveforms, t_target, n_pre::Integer, n_post::Integer)
    t = _shared_time_axis(wvfs)
    i_target = _sample_index(t, t_target)
    idxs = (i_target - n_pre):(i_target + n_post)
    firstindex(t) <= first(idxs) && last(idxs) <= lastindex(t) || throw(ArgumentError("window $idxs leaves the record $(axes(t, 1))"))
    stack(view(s, idxs) for s in wvfs.signal)
end
export aligned_windows

"""
    dplms_noise_matrix(wvfs::ArrayOfRDWaveforms, bl_window::ClosedInterval, n::Integer)

Noise matrix for a filter of `n` samples: the mean of `b bᵀ` over every window `b` of `n`
consecutive samples inside `bl_window` of every waveform. The waveforms must be baseline
subtracted; a residual offset appears as a rank-one term, not as noise. When the baseline
mean subtracted from the waveforms was estimated over the same window, the lowest-frequency
noise power is underestimated by about `n` divided by the number of baseline samples.
"""
function dplms_noise_matrix(wvfs::ArrayOfRDWaveforms, bl_window::ClosedInterval, n::Integer)
    t = _shared_time_axis(wvfs)
    i_from = _sample_index(t, leftendpoint(bl_window))
    i_to = _sample_index(t, rightendpoint(bl_window))
    n_win = i_to - i_from + 1 - n + 1
    n_win >= 1 || throw(ArgumentError("baseline window holds $(i_to - i_from + 1) samples, fewer than the filter length $n"))
    T = float(eltype(first(wvfs.signal)))
    N = zeros(T, n, n)
    W = Matrix{T}(undef, n_win, n)
    for s in wvfs.signal
        for j in 1:n, k in 1:n_win
            W[k, j] = s[i_from + k + j - 2]
        end
        mul!(N, W', W, true, true)
    end
    N ./ (n_win * length(wvfs))
end
export dplms_noise_matrix

"""
    dplms_reference(windows::AbstractMatrix, plateau::AbstractUnitRange)

Superpulse of the aligned `windows` (one column per event): the row mean divided by its mean
over the `plateau` rows, so the reference has unit amplitude after the rise. Returns
`(ref, plateau)` with the divisor.
"""
function dplms_reference(windows::AbstractMatrix, plateau::AbstractUnitRange{<:Integer})
    Base.require_one_based_indexing(windows)
    r = vec(mean(windows, dims = 2))
    p = mean(view(r, plateau))
    p > 0 || throw(ArgumentError("superpulse plateau is $p, expected a positive amplitude"))
    (ref = r ./ p, plateau = p)
end
export dplms_reference

"""
    dplms_flattop_matrix(windows_wide::AbstractMatrix, ff::Integer)

Flat-top matrix for a filter of `n = size(windows_wide, 1) - 2ff` samples, from aligned pulse
shapes (one column each, e.g. events or superpulses) widened by `ff` samples on each side. It
is the mean over columns and over the `2ff` shifts of the outer product of the column's first
difference, so `xᵀ F x` is the mean squared slope of the filter output around the pick-off.
"""
function dplms_flattop_matrix(windows_wide::AbstractMatrix, ff::Integer)
    ff >= 1 || throw(ArgumentError("flat-top widening ff must be at least 1, got $ff"))
    Base.require_one_based_indexing(windows_wide)
    n = size(windows_wide, 1) - 2ff
    n >= 1 || throw(ArgumentError("windows of $(size(windows_wide, 1)) samples are too short for widening $ff"))
    D = diff(windows_wide, dims = 1)
    F_full = D * D' ./ size(D, 2)
    F = zeros(eltype(F_full), n, n)
    for s in 0:(2ff - 1)
        F .+= view(F_full, (1 + s):(n + s), (1 + s):(n + s))
    end
    F ./ (2ff)
end
export dplms_flattop_matrix

"""
    dplms_binned_reference(windows::AbstractMatrix, bins::AbstractVector{<:Integer}, plateau::Real)

Reference constraints per drift-time bin: the superpulse `r_b` of every bin `b` in `bins`,
each divided by the common `plateau`, combined into `mat = mean_b r_b r_bᵀ` and
`vec = mean_b r_b`. Bins are numbered `1:maximum(bins)` and every bin must hold an event.
Returns `(mat, vec, refs)` with `refs` the per-bin references stacked into one column per bin.
"""
function dplms_binned_reference(windows::AbstractMatrix, bins::AbstractVector{<:Integer}, plateau::Real)
    Base.require_one_based_indexing(windows, bins)
    axes(windows, 2) == axes(bins, 1) || throw(DimensionMismatch("windows $(axes(windows, 2)) and bins $(axes(bins, 1)) differ"))
    all(>=(1), bins) || throw(ArgumentError("bin indices must be positive, got minimum $(minimum(bins))"))
    refs = map(1:maximum(bins)) do b
        members = findall(==(b), bins)
        isempty(members) && throw(ArgumentError("drift-time bin $b holds no event"))
        vec(mean(view(windows, :, members), dims = 2)) ./ plateau
    end
    (mat = mean(r * r' for r in refs), vec = mean(refs), refs = stack(refs))
end
export dplms_binned_reference

"""
    dplms_synthesis(noise, ref, flattop, binned; nm, za, ft, ct = 0)

Solve the DPLMS system

    (nm·noise + ref refᵀ + za·𝟙 + ft·flattop + ct·binned.mat) x = ref + ct·binned.vec

for the filter coefficients `x` along the window and rescale them so that `x ⋅ ref == 1`.
`𝟙` is the all-ones matrix, so `za` penalizes the squared coefficient sum. Returns

- `coeffs`: the impulse response, `reverse(x)`, for `ConvolutionFilter(DirectConvolution(), coeffs)`;
  the filter output at the last sample of a window equals `x ⋅ window`,
- `noise_out`: `sqrt(xᵀ noise x)`, the predicted output noise in units of the reference amplitude,
- `zero_area`: `sum(x)`.

`nm`, `za`, `ft` and `ct` act only through their ratios: scaling all four by the same factor
gives the same filter after normalization.
"""
function dplms_synthesis(noise::AbstractMatrix, ref::AbstractVector, flattop::AbstractMatrix, binned; nm::Real, za::Real, ft::Real, ct::Real = 0)
    n = length(ref)
    size(noise) == size(flattop) == (n, n) || throw(DimensionMismatch("noise $(size(noise)), flat-top $(size(flattop)) and reference length $n differ"))
    M = nm .* noise .+ ref * ref' .+ za .+ ft .* flattop
    v = collect(float.(ref))
    if !iszero(ct)
        size(binned.mat) == (n, n) && length(binned.vec) == n || throw(DimensionMismatch("binned reference $(size(binned.mat)), $(length(binned.vec)) does not match length $n"))
        M .+= ct .* binned.mat
        v .+= ct .* binned.vec
    end
    x = try
        Symmetric(M) \ v
    catch err
        err isa Union{LinearAlgebra.SingularException, LinearAlgebra.ZeroPivotException} || rethrow()
        throw(ArgumentError("DPLMS system is singular for nm=$nm, za=$za, ft=$ft, ct=$ct"))
    end
    d = dot(x, ref)
    iszero(d) && throw(ArgumentError("DPLMS solution has zero response to the reference for nm=$nm, za=$za, ft=$ft, ct=$ct"))
    x ./= d
    all(isfinite, x) || throw(ArgumentError("DPLMS solution is not finite for nm=$nm, za=$za, ft=$ft, ct=$ct"))
    (coeffs = reverse(x), noise_out = sqrt(dot(x, noise, x)), zero_area = sum(x))
end
export dplms_synthesis
