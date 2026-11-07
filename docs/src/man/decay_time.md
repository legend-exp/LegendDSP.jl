# Decay Time

The charge-sensitive preamplifier of an HPGe detector converts the collected charge into a voltage step that relaxes back to the baseline with an exponential tail,

```math
V(t) = V_0 + A \, e^{-(t - t_s)/\tau} \quad (t > t_s),
```

with the decay time ``\tau`` of the preamplifier. Pole-zero (PZ) correction cancels this tail so that the waveform after the rise stays flat. All energy filters of the HPGe chains ([HPGe DSP Chains](@ref)) are applied to the PZ-corrected waveform, so they need a good value of ``\tau``. The package determines ``\tau`` per waveform from the slope of the logarithm of the tail, and the production dataflow combines the values of many waveforms into one decay time per detector.

## Tail Statistics

[`tailstats`](@ref) fits a straight line to the natural logarithm of the samples in a time window ``[t_1, t_2]``:

```math
\ln y(t) = a + b\,t, \qquad \tau = -1/b.
```

It returns a named tuple with three fields:

| Field | Meaning |
|:------|:--------|
| `mean` | Mean of ``\ln y`` in the window |
| `sigma` | Standard deviation of ``\ln y`` in the window |
| `τ` | Decay time ``-1/b`` in the unit of the time axis |

The window must lie in the tail, after the rise has finished and before the waveform ends. The waveform must be baseline-subtracted, because the logarithm is taken of the sample values themselves. If any sample in the window is zero or negative, the logarithm is undefined and `tailstats` returns zero for all three fields. The estimate is biased when the noise is not small compared to the signal in the window, because the logarithm of a noisy sample is not symmetric around the logarithm of the mean; waveforms of large amplitude give the best estimates.

```@example decay_time
using LegendDSP
using RadiationDetectorSignals, RadiationDetectorDSP
using Unitful, StableRNGs, Statistics, ArraysOfArrays, IntervalSets

rng = StableRNG(7)
τ_true = 450.0u"µs"

function synthetic_waveform(amplitude; n = 8192, dt = 16.0u"ns", t_rise = 50.0u"µs", rise = 1.0u"µs", σ = 3.0)
    t = range(0.0u"ns", step = dt, length = n)
    x = ustrip.(u"µs", t) .- ustrip(u"µs", t_rise)
    shape = map(x) do xi
        xi < 0 ? 0.0 : xi < ustrip(u"µs", rise) ? xi / ustrip(u"µs", rise) :
            exp(-(xi - ustrip(u"µs", rise)) / ustrip(u"µs", τ_true))
    end
    RDWaveform(t, 1000 .+ amplitude .* shape .+ σ .* randn(rng, n))
end

wvfs = ArrayOfRDWaveforms([synthetic_waveform(4000 + 50 * i) for i in 1:50])
bl_window = 0.0u"µs" .. 39.0u"µs"
tail_window = 70.0u"µs" .. 110.0u"µs"

# Baseline-subtracted first waveform
wvf = shift_waveform(wvfs[1], -signalstats(wvfs[1], leftendpoint(bl_window), rightendpoint(bl_window)).mean)
stats = tailstats(wvf, leftendpoint(tail_window), rightendpoint(tail_window))
(τ = uconvert(u"µs", stats.τ), mean_of_log = stats.mean, sigma_of_log = stats.sigma)
```

## Decay Times of a Waveform Set

[`dsp_decay_times`](@ref) does the baseline subtraction and the tail fit for a set of waveforms. It has two methods:

```julia
dsp_decay_times(wvfs, bl_window, tail_window)
dsp_decay_times(wvfs, config::DSPConfig)
```

The second method takes the two windows from `config.bl_window` and `config.tail_window`. Both windows are closed intervals of Unitful times. The function decodes encoded waveforms, subtracts the mean of the samples in `bl_window` from each waveform, calls [`tailstats`](@ref) in `tail_window`, and returns the decay times as a vector of quantities in µs, one per waveform. Waveforms for which `tailstats` returns zero (non-positive samples in the window) give `0.0 µs`.

```@example decay_time
taus = dsp_decay_times(wvfs, bl_window, tail_window)
(median = median(taus), minimum = minimum(taus), maximum = maximum(taus))
```

The distribution of decay times of a set of waveforms is narrow and well described by a single peak. The scatter is driven by noise and by the amplitude of the waveforms. The production dataflow selects waveforms from a high-energy gamma line after a basic quality selection (see [`dsp_qc_flt_optimization_compressed`](@ref) and [Use in Juleana.jl](@ref)), computes `dsp_decay_times` on the presummed waveforms, restricts the values to a plausible interval, fits a truncated Gaussian to their histogram, and stores the center of the fit as the decay time of the detector.

The decay time is independent of the amplitude scale. For compressed raw data (see [Waveforms and Configuration](@ref)) the presummed waveform is passed to `dsp_decay_times`: the baseline subtraction and the logarithmic slope are the same as for the full-rate waveform, only the time axis step differs.

## Use in the Pole-Zero Correction

The decay time enters all chains as the argument `τ`; [`dsp_icpc`](@ref) and [`dsp_icpc_compressed`](@ref) build `InvCRFilter(τ)` from it and apply it to the baseline-subtracted waveforms before any energy filter. The same holds for the filter optimization scans ([Filter Optimization](@ref)).

`InvCRFilter(τ)` is the inverse of a CR high-pass filter with time constant ``\tau``. With the correct ``\tau`` the exponential tail turns into a constant, so that the slope of the waveform in the tail window vanishes. With a ``\tau`` that is too large the correction is too weak and the tail keeps a negative slope; with a ``\tau`` that is too small the corrected waveform rises in the tail. The chains write the mean, standard deviation and slope of the corrected tail window to `tailmean`, `tailsigma` and `tailslope` so that the quality of the correction can be checked for every event.

```@example decay_time
τ_est = median(taus)

function tail_slope(τ)
    corrected = InvCRFilter(τ).(shift_waveform.(wvfs, -signalstats.(wvfs, leftendpoint(bl_window), rightendpoint(bl_window)).mean))
    median(signalstats.(corrected, leftendpoint(tail_window), rightendpoint(tail_window)).slope)
end

[(τ = τ, slope = tail_slope(τ)) for τ in (0.9 * τ_est, τ_est, 1.1 * τ_est)]
```

The slope is in ADC counts per time unit of the time axis (ns). It changes sign around the estimated decay time.
