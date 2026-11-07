# Waveform Parameters

The functions on this page extract single numbers from waveforms: start times, threshold-crossing times, amplitudes of the current pulse, charge integrals, and summary statistics of a time window.
The DSP chains described in [HPGe DSP Chains](@ref) call them in sequence, and they can be used on their own.
The functions whose names start with `get_` take an `ArrayOfRDWaveforms` and return one value per waveform.
The statistics functions at the end of the page take one waveform and are broadcast over an array of waveforms.
Times are `Unitful` quantities.

```@setup routines
using CairoMakie
CairoMakie.activate!(type = "svg")
```

## Synthetic Waveforms

The examples use charge pulses that resemble the preamplifier output of a detector: a baseline offset, a rise at a start time `t0`, and an exponential decay with a decay time `tau`.
The noise is Gaussian with unit standard deviation.

```@example routines
using LegendDSP
using RadiationDetectorDSP
using RadiationDetectorSignals
using ArraysOfArrays
using IntervalSets
using StableRNGs
using Unitful

function decay_pulse(rng, t; t0, amplitude, rise = 0.3, tau = 50.0, offset = 10.0, noise = 1.0)
    x = ustrip.(u"µs", t)
    signal = [xi < t0 ? 0.0 : amplitude * (1 - exp(-(xi - t0) / rise)) * exp(-(xi - t0) / tau) for xi in x]
    RDWaveform(t, offset .+ signal .+ noise .* randn(rng, length(t)))
end

rng = StableRNG(3)
t = range(0.0u"µs", step = 0.01u"µs", length = 4001)
wvfs = ArrayOfRDWaveforms([
    decay_pulse(rng, t; t0 = t0, amplitude = a) for (t0, a) in zip((10.0, 12.0, 14.0), (50.0, 100.0, 150.0))
])
length(wvfs)
```

The first steps of the chains are applied here as well: the baseline is determined in a window before the pulse and subtracted, and the exponential decay is removed with the pole-zero correction `InvCRFilter` of RadiationDetectorDSP.jl:

```@example routines
bl_window = 0u"µs" .. 5u"µs"
bl_stats = signalstats.(wvfs, leftendpoint(bl_window), rightendpoint(bl_window))
wvfs_bl = shift_waveform.(wvfs, -bl_stats.mean)
wvfs_pz = InvCRFilter(50u"µs").(wvfs_bl)

round.(bl_stats.mean, digits = 2)
```

## Start Time and Threshold Crossings

[`get_t0`](@ref)`(wvfs_pz, t0_threshold; flt_pars, mintot)` determines the start time t0 of each pulse.
It filters the waveforms with an asymmetric trapezoidal filter, `TrapezoidalChargeFilter(flt_pars...)`, and takes the first time at which the output stays above `t0_threshold` for at least `mintot`.
The keyword `flt_pars` is a vector of the three filter times (averaging time, gap time, and averaging time on the other side of the gap); its default is `[40u"ns", 100u"ns", 2000u"ns"]`.
The default of `mintot` is `1500u"ns"`.
The threshold is in units of the signal.
The result is a vector of times in µs, with the value 0 µs for a waveform without a crossing.

[`get_threshold`](@ref)`(wvfs, threshold; mintot)` does the same without the filter: it intersects the waveforms themselves with `threshold`, which is either one number or an array with one threshold per waveform.
The default of `mintot` is `1000u"ns"`.
The chains use it for the times at which the pulse reaches a fraction of its height.

```@example routines
t0 = get_t0(wvfs_pz, 5.0)
t_20 = get_threshold(wvfs_pz, 20.0)
t_frac = get_threshold(wvfs_pz, [25.0, 50.0, 75.0])
(t0, t_20, t_frac)
```

Both functions interpolate linearly between samples.

## Current Pulse Maximum

The maximum of the current pulse, A, is the numerator of the A/E parameter.
[`get_a`](@ref)`(wvfs, current_window)` and [`get_wvf_maximum`](@ref)`(signal, start, stop)` return the maximum of a waveform within a time window as a `NamedTuple` `(max, t)`.
The maximum sample and its two neighbors define a parabola, and the height and the position of the vertex of the parabola are returned.
If the maximum sample is at the boundary of the window, the value and time of that sample are returned.
`get_a` applies `get_wvf_maximum` to each waveform of an array with the same window; its input is the current pulse, for example the output of one of the filters in [Filters](@ref).

```@example routines
wvfs_current = MovingWindowMultiFilter(200u"ns").(DerivativeFilter(1).(wvfs_pz))
a = get_a(wvfs_current, 5u"µs" .. 30u"µs")
(a.max, a.t)
```

```@example routines
get_wvf_maximum(wvfs_current[2], 5u"µs", 30u"µs")
```

## Drift Charge and LQ

[`get_qdrift`](@ref)`(wvfs, t_start, Δt; pol_power = 3, sign_est_length = 100u"ns")` is the difference of two charge integrals.
The waveforms are integrated first (`IntegratorFilter(1)`), and the value of the integral at a time is evaluated by a `SignalEstimator(PolynomialDNI(pol_power, sign_est_length))`, a polynomial fit of order `pol_power` to the integral over a window of length `sign_est_length` around the time.
With `Δt` a range of times, `first(Δt)` and `last(Δt)` are offsets relative to `t_start` (one start time per waveform), and the result is

```
(I(t_start + last(Δt)) - I(t_start + first(Δt))) - (I(t_start + first(Δt)) - I(t_start))
```

where `I` is the integral: the charge in the second part of the interval minus the charge in the first part.
`Δt` must be a `StepRangeLen` of times, for example `range(1u"µs", 2u"µs", length = 2)`, and only its first and last element are used.

[`get_lq`](@ref)`(wvfs_int, t_start, early, gap, late; pol_power, sign_est_length)` generalizes this with separate windows.
It takes waveforms that are already integrated, so that the integration is done once for several LQ windows, and returns the charge in the late window `[t_start + early + gap, t_start + early + gap + late]` minus the charge in the early window `[t_start, t_start + early]`.
`early`, `gap`, and `late` are times.
`get_qdrift(wvfs, t, Δt)` is equal to `get_lq(wvfs_int, t, first(Δt), 0u"µs", last(Δt) - first(Δt))`.

```@example routines
wvfs_int = IntegratorFilter(1).(wvfs_pz)

q_drift = get_qdrift(wvfs_pz, t0, range(1u"µs", 2u"µs", length = 2))
lq = get_lq(wvfs_int, t0, 1u"µs", 0u"µs", 1u"µs")
(q_drift, lq)
```

## Triggers in the Waveform

[`get_triggers`](@ref)`(wvfs, sigma_threshold, bl_window; mintot = 100.0u"ns")` finds crossings of a threshold that is a multiple of the baseline noise.
The standard deviation of each waveform in `bl_window` (counted from the first time of the first waveform) is multiplied by `sigma_threshold`, and the waveform is intersected with this threshold.
If the standard deviation is zero, as for a saturated waveform, the threshold is 1.
The waveforms must have a baseline of about zero, because the threshold is not relative to the baseline mean.
The result is a `NamedTuple` `(intersect, n)` with the time of the first crossing and the number of crossings in each waveform; `intersect` is 0 µs if there is no crossing.
The number of crossings is used to detect pile-up inside a trace.

[`get_triggers_reversed`](@ref) does the same on the time-reversed waveform and returns the time in the original orientation, `last(time) - x`.
It therefore locates the last time at which the waveform falls below the threshold; the result is the last time of the waveform if the waveform does not end below the threshold.
The pulses above decay slowly and are still above the threshold at the end of the trace, so the example for this function uses a pulse with a decay time of 2 µs.
`n` also counts crossings that noise causes where the pulse falls to the threshold.

```@example routines
trig = get_triggers(wvfs_bl, 5, bl_window)
(trig.intersect, trig.n)
```

```@example routines
wvfs_fast = ArrayOfRDWaveforms([decay_pulse(rng, t; t0 = 10.0, amplitude = 100.0, tau = 2.0, offset = 0.0)])
trig_fwd = get_triggers(wvfs_fast, 5, bl_window)
trig_rev = get_triggers_reversed(wvfs_fast, 5, bl_window)
(trig_fwd.intersect, trig_rev.intersect, trig_rev.n)
```

## Extrema, Tail, Threshold, and Saturation Statistics

The functions in this section take a single waveform, either an `RDWaveform` with a time window given as times, or a plain vector of samples (`AbstractSamples`) with a window given as sample indices.
Each returns a `NamedTuple`, and broadcasting over an array of waveforms returns the fields as arrays.

[`extremestats`](@ref)`(signal, start, stop)` returns the minimum and maximum sample and their times, `(min, max, tmin, tmax)`.
Without `start` and `stop` it uses the whole signal.
The values are those of the samples; no interpolation is done.

```@example routines
es = extremestats.(wvfs_bl)
(es.min, es.max, es.tmax)
```

[`tailstats`](@ref)`(signal, start, stop)` fits a straight line to the natural logarithm of the samples in the window and returns `(mean, sigma, τ)`: the mean and the standard deviation of the logarithm of the samples and the decay time `τ = -1/slope` of the fit.
The signal must be positive in the whole window, otherwise all three values are zero.
For a waveform without baseline offset the decay time is that of the exponential tail.
It is the basis of the decay-time determination in [Decay Time](@ref).

```@example routines
tail = tailstats.(wvfs_bl, 20u"µs", 40u"µs")
tail.τ
```

[`thresholdstats`](@ref)`(signal, min, max)` returns the standard deviation of those samples of `signal` that lie between `min` and `max`.
The bounds are optional; a missing bound is ignored.
It estimates the baseline noise when the bounds exclude the pulse.
[`thresholdstats_mad`](@ref) uses the median absolute deviation of the same samples, multiplied by 1.4826 so that it equals the standard deviation of normally distributed data, which is less sensitive to a few outliers.
It returns 0 if no sample lies within the bounds, whereas `thresholdstats` returns `NaN` in that case.

```@example routines
(thresholdstats(wvfs_bl[1], -5.0, 15.0), thresholdstats_mad(wvfs_bl[1], -5.0, 15.0), thresholdstats(wvfs_bl[1]))
```

[`saturation`](@ref)`(signal, low, high)` and `saturation(signal, start, stop, low, high)` count the samples that are equal to the digitizer limits `low` and `high`.
It returns `(low, high, max_cons_low, max_cons_high)`: the number of samples at each limit and the longest run of consecutive samples at each limit.
A pulse that exceeds the input range of the digitizer stays at `high` for several samples.

```@example routines
clipped = RDWaveform(t, clamp.(wvfs[3].signal, 0.0, 100.0))
saturation(clipped, 0.0, 100.0)
```

```@example routines
using CairoMakie

fig = Figure(size = (700, 300))
ax = Axis(fig[1, 1], xlabel = "Time (µs)", ylabel = "Signal (arb. units)")
lines!(ax, ustrip.(u"µs", wvfs[3].time), wvfs[3].signal, label = "waveform")
lines!(ax, ustrip.(u"µs", clipped.time), clipped.signal, label = "clipped at 100")
vlines!(ax, ustrip(u"µs", t0[3]), color = :gray, linestyle = :dash, label = "get_t0")
axislegend(ax)
fig
```
