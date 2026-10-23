# Filters

LegendDSP.jl adds a set of filter types to those of [RadiationDetectorDSP.jl](https://github.com/JuliaPhysics/RadiationDetectorDSP.jl).
They are used to estimate the current pulse from a charge pulse (needed for the A/E parameter), to smooth waveforms, to reduce the sampling rate, to change the time axis, and to find the times at which a signal crosses a threshold.
This page describes how the filters are constructed and applied and what each of them does.
The functions that call them to produce DSP parameters are described in [HPGe DSP Chains](@ref) and [Waveform Parameters](@ref).

```@setup filters
using CairoMakie
CairoMakie.activate!(type = "svg")
```

## Constructing and Applying Filters

A filter is a plain immutable object that holds the filter parameters.
Parameters with the dimension of a time or a frequency are Unitful quantities, so that the same filter can be applied to waveforms with different sampling periods:

```@example filters
using LegendDSP
using RadiationDetectorSignals
using Unitful
using StableRNGs

flt_mwa = MovingWindowMultiFilter(200u"ns"; passes = 3)
```

A filter is applied by calling it on a waveform, [`RDWaveform`](https://juliaphysics.github.io/RadiationDetectorSignals.jl/stable/), which returns a new waveform.
To process many waveforms, broadcast the filter over an `ArrayOfRDWaveforms`; the result is again an `ArrayOfRDWaveforms`.
The filter parameters in units of time are converted to a number of samples using the step of the time axis of the waveform, so the time axis must be an `AbstractRange`.

The examples on this page use a synthetic charge pulse: a logistic step of height 100 centered at 8 µs with a width parameter of 0.3 µs, sampled every 10 ns and superimposed with Gaussian noise of unit standard deviation.
A second set of waveforms with a different noise realization shows the broadcasting form.

```@example filters
function synthetic_pulse(rng, t; t0 = 8.0, amplitude = 100.0, rise = 0.3, noise = 1.0)
    signal = amplitude ./ (1 .+ exp.(-(ustrip.(u"µs", t) .- t0) ./ rise)) .+ noise .* randn(rng, length(t))
    RDWaveform(t, signal)
end

rng = StableRNG(1)
t = range(0.0u"µs", step = 0.01u"µs", length = 2001)
wvf = synthetic_pulse(rng, t)
wvfs = ArrayOfRDWaveforms([synthetic_pulse(rng, t) for _ in 1:4])

wvf_smooth  = flt_mwa(wvf)
wvfs_smooth = flt_mwa.(wvfs)
typeof(wvfs_smooth) == typeof(wvfs), length(wvf_smooth.signal)
```

Filters that act sample by sample on the signal (for example [`DerivativeFilter`](@ref) and [`MovingWindowFilter`](@ref)) return a waveform with the same time axis.
Filters that change the number of samples or the time axis ([`HaarAveragingFilter`](@ref), [`TimeAxisFilter`](@ref)) return a waveform with a new time axis.
The filter types of RadiationDetectorDSP.jl are applied in the same way.
Their interface, which consists of `fltinstance` and `rdfilt!` for the low-level implementation, is described in the RadiationDetectorDSP.jl documentation.

## Current Estimation

The current pulse of a detector is the time derivative of the charge pulse.
The A/E parameter uses the maximum of the current pulse, A, and therefore needs a derivative that is smoothed enough to suppress the noise of the sampled charge pulse.
LegendDSP.jl offers two ways to obtain it:

- `SavitzkyGolayFilter(length, degree, 1)` from RadiationDetectorDSP.jl fits a polynomial of the given `degree` in a sliding window of the given `length` and returns its first derivative.
- [`DerivativeFilter`](@ref) followed by [`MovingWindowMultiFilter`](@ref) takes the difference of consecutive samples and smooths it with several moving average passes.

```@example filters
using RadiationDetectorDSP

wvf_raw_deriv = DerivativeFilter(1)(wvf)
wvf_mwa_deriv = MovingWindowMultiFilter(400u"ns"; passes = 3)(wvf_raw_deriv)
wvf_sg_deriv  = SavitzkyGolayFilter(1u"µs", 3, 1)(wvf)

extrema(wvf_raw_deriv.signal), maximum(wvf_mwa_deriv.signal), maximum(wvf_sg_deriv.signal)
```

[`DerivativeFilter`](@ref) computes `gain * (x[i] - x[i-1])` and does not divide by the sampling period, so its output is in units of the signal per sample.
The first sample uses `x[2] - x[1]`.
The maximum of the derivative of the synthetic pulse is `100 / (4 * 0.3 µs)`, which is 0.83 per 10 ns sample.
The `SavitzkyGolayFilter` derivative has the unit of the signal per sample as well, and it drops the samples at both ends of the waveform for which the window does not fit; its output waveform is therefore shorter than the input.

```@example filters
function plot_waveforms(pairs...; xlims, ylabel)
    fig = Figure(size = (700, 320))
    ax = Axis(fig[1, 1], xlabel = "Time (µs)", ylabel = ylabel, limits = (xlims, nothing))
    for (label, w) in pairs
        lines!(ax, ustrip.(u"µs", w.time), w.signal, label = label)
    end
    axislegend(ax, position = :lt)
    fig
end

plot_waveforms(
    "DerivativeFilter" => wvf_raw_deriv,
    "+ MovingWindowMultiFilter" => wvf_mwa_deriv,
    "SavitzkyGolayFilter" => wvf_sg_deriv;
    xlims = (6, 11), ylabel = "Current (arb. units)",
)
```

In the DSP chains, the window length of either estimator is the free parameter that is optimized with a scan over a grid of lengths (see [Filter Optimization](@ref)).

## Smoothing Filters

The following filters smooth the signal itself.
They have the same number of output and input samples.

[`MovingWindowFilter`](@ref) is a causal moving average over `length`: each output sample is the mean of the input samples in the window that ends at that sample, so the output lags the input by half a window.
For the first `length` samples the window is padded with the first input sample.
[`MovingWindowMultiFilter`](@ref) applies the same average `passes` times (default 3) and alternates between forward and backward passes.
An even number of passes has no net delay; with an odd number of passes the remaining delay is that of a single pass.
`passes` must be at least 1, otherwise an `ArgumentError` is thrown.

[`WeightedSavitzkyGolayFilter`](@ref) is a Savitzky-Golay smoother of polynomial `degree` whose window is weighted by a function selected with `weightType`:

| `weightType` | Weight function |
|:---|:---|
| 0 | none (unweighted Savitzky-Golay filter) |
| 1 | Gaussian-like |
| 2 | Hann |
| 3 | squared Hann |
| 4 | cubed Hann |

Other values of `weightType` are rejected by the constructor.
The `length` is the full window length.

[`ModifiedSincFilter`](@ref) is a finite impulse response low-pass filter with a modified sinc kernel.
`d` is an even degree between 2 and 10, which sets the number of extrema in the kernel, and `m` is the half-width of the kernel as a time.
The kernel must contain at least `d / 2 + 2` samples on each side; otherwise an error is thrown when the filter is applied.
The default `m = 3` has no unit and therefore cannot be used with a time axis with units; pass a time explicitly.

[`WhittakerHendersonFilter`](@ref) returns the signal that minimizes the sum of the squared deviations from the input plus `λ` times the sum of the squared differences of order `p`.
Both parameters are dimensionless and refer to differences between neighboring samples, so the same `λ` smooths more or less depending on the sampling period of the waveform.
A larger `λ` gives a smoother output.

```@example filters
flt_wsg = WeightedSavitzkyGolayFilter(length = 500u"ns", degree = 2, weightType = 2)
flt_sinc = ModifiedSincFilter(d = 4, m = 1u"µs")
flt_wh = WhittakerHendersonFilter(p = 3, λ = 1000)

plot_waveforms(
    "raw" => wvf,
    "MovingWindowMultiFilter" => flt_mwa(wvf),
    "WeightedSavitzkyGolayFilter" => flt_wsg(wvf),
    "ModifiedSincFilter" => flt_sinc(wvf),
    "WhittakerHendersonFilter" => flt_wh(wvf);
    xlims = (7, 9), ylabel = "Charge (arb. units)",
)
```

[`WeightedSavitzkyGolayFilter`](@ref) is used for the smoothing of PMT waveforms in [`dsp_pmts`](@ref) (see [SiPM and PMT Processing](@ref)).
The other smoothing filters are available for studies of the filter response and are not called by the DSP chains.

## Changing the Sampling Rate and the Time Axis

[`HaarAveragingFilter`](@ref)`(down_sampling_rate)` sums two consecutive samples, scales the sum by `1/sqrt(2)`, and keeps one result per `down_sampling_rate` input samples.
With `down_sampling_rate = 2` each application halves the number of samples; a constant signal level is multiplied by `sqrt(2)` per application.
The time axis of the output is `time[1:down_sampling_rate:end]`.
The classifiers for the machine-learning quality cut apply this filter repeatedly (see [Machine-Learning Quality Cuts](@ref)).

[`TimeAxisFilter`](@ref)`(period, offset)` leaves the signal unchanged and replaces the time axis by a range that starts at `first(time) + offset` and has step `period`.
[`dsp_pmts`](@ref) uses it to set the sampling period of PMT waveforms to the value given in its configuration.
`offset` defaults to zero, and both arguments must have the same type.
A time axis that is not a range raises an `ArgumentError`.

```@example filters
wvf_haar = HaarAveragingFilter(2)(wvf)
wvf_axis = TimeAxisFilter(4u"ns")(wvf)

length(wvf.signal), length(wvf_haar.signal), step(wvf_haar.time),
step(wvf_axis.time), length(wvf_axis.signal)
```

## Finding Threshold Crossings

Two callable types find the times at which a signal crosses a threshold.
They use the `Intersect` filter of RadiationDetectorDSP.jl as the model, which returns the first crossing, and extend it to several crossings per waveform and to several thresholds.

[`IntersectMaximum`](@ref) finds all up-crossings of a threshold that stay above it for at least `mintot` and searches for the maximum of the signal in the window of length `maxtot` that follows each crossing.
It is called as `IntersectMaximum(mintot, maxtot)(waveform, threshold)` and returns a `NamedTuple`:

| Field | Content |
|:---|:---|
| `x` | times of the up-crossings, linearly interpolated between samples |
| `x_high` | times of the following down-crossings; the last time of the waveform if the signal stays above the threshold |
| `x_tot` | `x_high - x`, the time over threshold |
| `max` | maxima in the windows, interpolated with a parabola through the three samples around the maximum |
| `multiplicity` | number of crossings |

`mintot` and `maxtot` have the type of the elements of the time axis, so both must have the same type: `IntersectMaximum(100u"ns", 1000u"ns")` works, while `IntersectMaximum(100u"ns", 1u"µs")` does not, because the units differ.
The default values, 4 and 100, have no unit and are suitable only for waveforms with a dimensionless time axis.
The type is used for trigger finding on SiPM and PMT waveforms (see [SiPM and PMT Processing](@ref)).
The following example finds two pulses in a synthetic current waveform.

```@example filters
rng = StableRNG(2)
wvf_pileup = let w1 = synthetic_pulse(rng, t; t0 = 5.0, amplitude = 100.0),
                 w2 = synthetic_pulse(rng, t; t0 = 14.0, amplitude = 40.0, noise = 0.0)
    RDWaveform(t, w1.signal .+ w2.signal)
end

current = MovingWindowMultiFilter(200u"ns")(DerivativeFilter(1)(wvf_pileup))
threshold = 8 * thresholdstats_mad(current)

triggers = IntersectMaximum(100u"ns", 1000u"ns")(current, threshold)
triggers.multiplicity, round.(u"µs", triggers.x, digits = 2), round.(triggers.max, digits = 2)
```

`thresholdstats_mad` is described in [Waveform Parameters](@ref).

[`MultiIntersect`](@ref) returns for every entry of `threshold_ratios` the time at which the signal first exceeds that fraction of its maximum, for a time over threshold of at least `mintot`.
Each crossing is refined by a polynomial of degree `d` that is fitted to the `2n` samples around it, evaluated on a grid that is upsampled by `sampling_rate`.
The result is a vector with one time per ratio.
The default ratios are `0.01:0.01:0.9`.
The time axis of the waveform must have a floating-point or integer element type; an axis with rational elements, as produced by `range(0u"µs", 20u"µs", step = 10u"ns")`, raises an error.
A threshold for which no crossing is found in the fit window gives the value 0, and an `AssertionError` is thrown if a crossing is closer than `n` samples to the start or the end of the waveform.

```@example filters
flt_multi = MultiIntersect(threshold_ratios = [0.1, 0.5, 0.9], mintot = 100u"ns")
t10, t50, t90 = flt_multi(wvf)
round.(u"µs", (t10, t50, t90), digits = 2)
```

## Filters of RadiationDetectorDSP.jl Used by the DSP Chains

The DSP chains also use the following filter types of RadiationDetectorDSP.jl.
They are not defined in LegendDSP.jl; see the documentation of RadiationDetectorDSP.jl for their parameters.

| Filter | Use in LegendDSP.jl |
|:---|:---|
| `InvCRFilter(τ)` | pole-zero correction of the exponential decay of the preamplifier output with decay time `τ` |
| `TrapezoidalChargeFilter(rise, flat_top)` | energy estimate with a trapezoidal filter; also the fast asymmetric trapezoid in [`get_t0`](@ref) |
| `CUSPChargeFilter` | energy estimate with a cusp-shaped weighting function |
| `ZACChargeFilter` | energy estimate with a zero-area cusp-shaped weighting function |
| `SavitzkyGolayFilter(length, degree, derivative)` | smoothing and derivative of the waveform |
| `TruncateFilter(interval)` | restriction of a waveform to a time interval |
| `IntegratorFilter` | integration of a current or derivative waveform |
| `Intersect(mintot = ...)` | first crossing of a threshold, used by [`get_t0`](@ref), [`get_threshold`](@ref) and [`get_triggers`](@ref) |
| `SignalEstimator(PolynomialDNI(order, length))` | value of a filtered waveform at a given time, interpolated with a polynomial |

```@example filters
flt_trap = TrapezoidalChargeFilter(1000u"ns", 2000u"ns")
wvf_trap = flt_trap(wvf)
maximum(wvf_trap.signal)
```

The flat top of the trapezoid reaches the height of the step, here 100, when it is longer than the rise of the pulse.
