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
