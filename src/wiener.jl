# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    wiener_template(windows::AbstractMatrix, amplitudes::AbstractVector)

Template of the aligned `windows` (one column per event): every window is divided by
its amplitude, the median over the events is taken per row, and the result is divided by its
maximum, so the template has the height 1. The amplitudes need a common scale only: any
quantity proportional to the pulse height serves. They should not be the maxima of the
windows themselves, which the noise raises.

Throws if the numbers of windows and amplitudes differ, if an amplitude is not positive or
if the maximum of the median is not positive.

# Returns
- named tuple `(template, peak)` with the template and the divisor, the pulse height per unit of amplitude
"""
function wiener_template(windows::AbstractMatrix, amplitudes::AbstractVector)
    Base.require_one_based_indexing(windows, amplitudes)
    axes(windows, 2) == axes(amplitudes, 1) || throw(DimensionMismatch("windows $(axes(windows, 2)) and amplitudes $(axes(amplitudes, 1)) differ"))
    all(a -> a > zero(a), amplitudes) || throw(ArgumentError("amplitudes must be positive, got minimum $(minimum(amplitudes))"))
    r = vec(median(windows ./ amplitudes', dims = 2))
    p = maximum(r)
    p > zero(p) || throw(ArgumentError("template maximum is $p, expected a positive amplitude"))
    template = r ./ p
    (template = template, peak = p)
end
export wiener_template

"""
    wiener_noise_spectrum(traces::AbstractMatrix{<:Real}, dt::RealQuantity)

One-sided noise power spectrum of the `traces` (one column per trace, sampled every `dt`):
the median over the traces of their Hann-windowed periodograms, divided by `log(2)`, which
turns the median of an exponentially distributed periodogram value into its mean. A pulse in
a few of the traces hardly moves the median. `sum(psd) * step(freq)` is the variance of the
samples. The factor `log(2)` holds for many traces: with `m` traces the spectrum comes out
high by about `0.72 / m`.

The traces must be baseline subtracted and have an even number of samples `n`, so that
`n = 2 (length(psd) - 1)`. The first and the last value, at zero frequency and at the Nyquist
frequency, are not noise estimates like the others: the Fourier transform is real there, so
the division by `log(2)` leaves them about a third low, and the baseline subtraction removes
power at the lowest frequencies. The expected squared magnitude of the discrete Fourier
transform (`rfft`) of a trace of `n` samples is `psd * n / (2 dt)` for the values in between
and `psd * n / dt` for the first and the last.

Throws if the number of samples is odd.

# Returns
- named tuple `(freq, psd)`: the frequencies `(0:n÷2) / (n dt)` and the power per frequency, in the units of the squared samples times the unit of `dt`
"""
function wiener_noise_spectrum(traces::AbstractMatrix{<:Real}, dt::RealQuantity)
    Base.require_one_based_indexing(traces)
    n = size(traces, 1)
    iseven(n) || throw(ArgumentError("noise traces have $n samples, expected an even number"))
    T = float(eltype(traces))
    w = [sinpi(T(k) / (n - 1))^2 for k in 0:n-1]
    P = abs2.(rfft(traces .* w, 1))
    P[2:end-1, :] .*= 2
    freq = (0:n÷2) ./ (n * dt)
    psd = vec(median(P, dims = 2)) .* (dt / (sum(abs2, w) * log(T(2))))
    (freq = freq, psd = psd)
end
export wiener_noise_spectrum

"""
    wiener_synthesis(noise, template::AbstractVector{<:Real}, target::AbstractVector{<:Real}, n_pre::Integer; λ::Real, n_half::Integer)

Wiener deconvolution kernel that turns the pulse shape `template` into `target`:

    H = G P* / (|P|² + λ N)

with `P` and `G` the Fourier transforms of `template` and `target` and `N` the noise power, all on
the grid of the noise traces. Where the pulse is strong the kernel divides its shape out and
puts the target in; where the noise is stronger than `|P|² / λ` the kernel goes to zero. `λ`
weighs the noise against a pulse of the height of `template`: a larger `λ` gives less noise after
the kernel and a response that follows the target less closely.

`template` and `target` share one axis with `n_pre` samples before the sample of lag 0, and are
taken as zero outside it, so a pulse that has not decayed at the end of `template` must be continued
by the caller. The kernel is cut to `n_half` samples on either side of lag 0, and `response`
and `noise_out` are those of the cut kernel. With a `target` of sum 1 the filtered trace has,
for every pulse, the target scaled to an area equal to the pulse height in units of `template`.

Throws if `template` and `target` differ in length or are longer than the noise traces, if `n_pre`
or `n_half` do not fit, or if the kernel is not finite.

# Arguments
- `noise`: noise spectrum from [`wiener_noise_spectrum`](@ref)
- `template`: pulse shape, for example from [`wiener_template`](@ref)
- `target`: the shape the kernel makes of `template`, on the axis of `template`
- `n_pre`: number of samples of `template` and `target` before lag 0

# Keywords
- `λ`: weight of the noise
- `n_half`: number of kernel samples on either side of lag 0

# Returns
- named tuple with

- `coeffs`: the `2 n_half + 1` kernel values, lag 0 in the middle; sample `i` of the filtered
  trace is the sum of `coeffs[n_half + 1 + l]` times sample `i - l` of the trace,
- `response`: the filtered `template`, on the axis of `template`,
- `noise_out`: the predicted standard deviation of a sample of the filtered noise, in the units of the noise traces,
- `spectrum`: named tuple of what the kernel is built from, on the frequencies `freq` of the noise spectrum: the squared
  magnitudes of the Fourier transforms of `template` and `target`, the expected one of a noise trace (`noise`, the `N` above)
  and the `gain` of the cut kernel, `|H|`.
"""
function wiener_synthesis(noise, template::AbstractVector{<:Real}, target::AbstractVector{<:Real}, n_pre::Integer; λ::Real, n_half::Integer)
    Base.require_one_based_indexing(template, target, noise.psd)
    m, n = length(template), 2 * (length(noise.psd) - 1)
    length(target) == m || throw(DimensionMismatch("template length $m and target length $(length(target)) differ"))
    0 <= n_pre < m <= n || throw(ArgumentError("template of $m samples with $n_pre before lag 0 does not fit the $n samples of the noise traces"))
    0 <= n_half < n ÷ 2 || throw(ArgumentError("kernel half length $n_half does not fit the $n samples of the noise traces"))
    on_grid(y) = circshift(vcat(float.(y), zeros(float(eltype(y)), n - m)), -n_pre)
    P, G = rfft(on_grid(template)), rfft(on_grid(target))
    df = step(noise.freq)
    N = ustrip.(NoUnits, noise.psd .* (df * n^2 / 2))
    N[1] *= 2
    N[end] *= 2
    h = irfft(G .* conj.(P) ./ (abs2.(P) .+ λ .* N), n)
    all(isfinite, h) || throw(ArgumentError("Wiener kernel is not finite for λ=$λ"))
    h[n_half+2:n-n_half] .= 0
    H = rfft(h)
    coeffs = circshift(h, n_half)[1:2n_half+1]
    response = circshift(irfft(H .* P, n), n_pre)[1:m]
    noise_out = sqrt(ustrip(NoUnits, sum(abs2.(H) .* noise.psd) * df))
    spectrum = (freq = noise.freq, template = abs2.(P), target = abs2.(G), noise = N, gain = abs.(H))
    (coeffs = coeffs, response = response, noise_out = noise_out, spectrum = spectrum)
end
export wiener_synthesis
