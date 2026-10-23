# SiPM and PMT Processing

The silicon photomultiplier (SiPM) and photomultiplier tube (PMT) channels of the liquid-argon instrumentation record short pulses of single photoelectrons on top of a baseline.
LegendDSP.jl provides one DSP chain for each of the two channel types, and a set of functions that scan the free parameters of the SiPM chain.
The chains of the HPGe detectors take a typed configuration object, [`DSPConfig`](@ref), which is described in [Waveforms and Configuration](@ref).
The SiPM and PMT functions on this page take the configuration as a plain `PropDict` of the PropDicts.jl package.
A `PropDict` holds nested key-value pairs; the keys that each function reads are listed below, and a missing key raises an error when the function is called.
Values with a unit, such as times, are stored as `Unitful` quantities in the dictionary.

```@setup sipm_pmt
using CairoMakie
CairoMakie.activate!(type = "svg")
```

## Input Tables and Synthetic Data

All chains take a `TypedTables.Table` with one row per event, call `decode_data` on it, and return a new `Table` with one row per event.
The examples below use a table of synthetic SiPM-like waveforms: a baseline of 20 ADC counts, Gaussian noise, and pulses with an instantaneous rise and an exponential decay.

```@example sipm_pmt
using LegendDSP
using ArraysOfArrays
using PropDicts
using RadiationDetectorSignals
using StableRNGs
using TypedTables
using Unitful

function sipm_waveform(rng, t; pulses, noise = 0.5, tau = 0.3)
    x = ustrip.(u"µs", t)
    signal = zeros(length(x))
    for (t0, amplitude) in pulses
        signal .+= [xi < t0 ? 0.0 : amplitude * exp(-(xi - t0) / tau) for xi in x]
    end
    RDWaveform(t, 20.0 .+ signal .+ noise .* randn(rng, length(x)))
end

rng = StableRNG(5)
t = range(0.0u"µs", step = 0.016u"µs", length = 2000)
pulses = [[(5.0, 10.0)], [(5.0, 10.0), (15.0, 20.0)], Tuple{Float64,Float64}[]]
wvfs = ArrayOfRDWaveforms([sipm_waveform(rng, t; pulses = p) for p in pulses])

data = Table(
    waveform = wvfs, waveform_bit_drop = wvfs, baseline = [20, 20, 20],
    timestamp = [1.0, 2.0, 3.0], eventnumber = [1, 2, 3], daqenergy = [100, 200, 300], channel = [1, 1, 1],
)
length(data)
```

## SiPM Chain

[`dsp_sipm`](@ref)`(data, config, pars_optimization)` and [`dsp_sipm_compressed`](@ref)`(data, config, pars_optimization)` process the waveforms of one SiPM channel.
The two functions have the same steps.
They differ only in the input column that holds the waveform: `waveform` for `dsp_sipm` and `waveform_bit_drop` for `dsp_sipm_compressed`.
The other input columns are `baseline`, `timestamp`, `eventnumber`, and `daqenergy`.

The steps are:

1. The waveforms are converted to floating-point values and the maximum and minimum with their times are determined for the whole trace and for the part of the trace within `t0_hpge_window`, the time interval in which the liquid-argon light is expected.
2. A Savitzky-Golay filter of window length `pars_optimization.sg.wl`, polynomial degree `sg_flt_degree` and derivative order 1 yields the smoothed derivative of the waveform. This is the SG branch.
3. The SG threshold is `n_σ_threshold` times the median-absolute-deviation estimate of the noise ([`thresholdstats_mad`](@ref)) of the samples of the derivative that lie between `min_threshold` and `max_threshold`. [`IntersectMaximum`](@ref) finds all up-crossings of this threshold and the maximum of the derivative after each crossing.
4. The derivative is integrated again (`IntegratorFilter`), which gives a baseline-corrected waveform. The baseline statistics are determined on the integrated waveform from its start to the first trigger, but at least over three samples. The waveform statistics are determined over the whole integrated waveform (`signalstats`).
5. Discharges are detected on the integrated waveform multiplied by -1. The threshold is `n_σ_dc_threshold` times the median-absolute-deviation estimate of the samples between `min_dc_threshold` and `max_dc_threshold`, and the crossings are found with [`IntersectMaximum`](@ref).
6. The trap branch applies the pole-zero correction `InvCRFilter(pz_tau)` and a `TrapezoidalChargeFilter(rt, ft)` to the integrated waveform. Triggers are found as in the SG branch with the thresholds of the trap section of the configuration.
7. The discharge detection of the trap branch is made on the flipped integrated waveform with its own thresholds. It uses the [`IntersectMaximum`](@ref) object that is built from the settings of the SG section (`min_tot_intersect` and `max_tot_intersect`).

### Configuration

`config` must provide the keys in the following list; the bullets show `PropDict` paths.

- `sg_flt_degree`: polynomial degree of the Savitzky-Golay filter.
- `t0_hpge_window`: pair of times (start, stop) for the liquid-argon window.
- `filters.sg.min_tot_intersect`, `filters.sg.max_tot_intersect`: minimum time over threshold and the window in which the maximum is searched, as times of the same unit.
- `filters.sg.min_threshold`, `filters.sg.max_threshold`: bounds of the samples that enter the noise estimate of the SG derivative.
- `filters.sg.n_σ_threshold`: factor between noise estimate and threshold.
- `filters.sg.min_dc_threshold`, `filters.sg.max_dc_threshold`, `filters.sg.n_σ_dc_threshold`: the same three values for the discharge detection.
- `filters.trap.rt`, `filters.trap.ft`: rise time and flat-top time of the trapezoidal filter.
- `filters.trap.pz_tau`: decay time of the pole-zero correction.
- `filters.trap.min_tot_intersect`, `filters.trap.max_tot_intersect`, `filters.trap.min_threshold`, `filters.trap.max_threshold`, `filters.trap.n_σ_threshold`, `filters.trap.min_dc_threshold`, `filters.trap.max_dc_threshold`, `filters.trap.n_σ_dc_threshold`: the corresponding values for the trap branch.

`pars_optimization` provides the optimized window length as `pars_optimization.sg.wl`.
The scan that determines it is described [below](@ref "SiPM Parameter Scans").

### Output

The output table has the following columns.

| Column | Content |
|:---|:---|
| `blfc`, `timestamp`, `eventID_fadc`, `e_fc` | the input columns `baseline`, `timestamp`, `eventnumber`, and `daqenergy` |
| `t_max`, `t_min` | times of the maximum and minimum of the waveform, in µs |
| `t_max_lar`, `t_min_lar` | the same within `t0_hpge_window`, in µs |
| `e_max`, `e_min`, `e_max_lar`, `e_min_lar` | maximum and minimum amplitudes, for the whole trace and in the liquid-argon window |
| `blmean`, `blsigma`, `blslope`, `bloffset` | baseline statistics of the integrated waveform |
| `wfmean`, `wfsigma`, `wfslope`, `wfoffset` | statistics of the whole integrated waveform |
| `threshold_sg`, `threshold_DC_sg` | thresholds of the SG branch for triggers and discharges |
| `trig_pos_sg`, `trig_max_sg` | positions and maxima of the SG triggers (vectors per event) |
| `trig_pos_DC_sg`, `trig_max_DC_sg` | positions and maxima of the SG discharge triggers |
| `threshold_trap`, `threshold_DC_trap` | thresholds of the trap branch |
| `trig_pos_trap`, `trig_pos_high_trap`, `trig_pos_tot_trap`, `trig_max_trap` | position, down-crossing position, time over threshold, and maximum of the trap triggers |
| `trig_pos_DC_trap`, `trig_pos_high_DC_trap`, `trig_pos_tot_DC_trap`, `trig_max_DC_trap` | the same for the trap discharge triggers |

The trigger columns hold one vector per event, stored as a `VectorOfVectors`; the number of triggers is the length of the vector.

```@example sipm_pmt
sipm_cfg = PropDict(
    :sg_flt_degree => 2,
    :t0_hpge_window => (10u"µs", 20u"µs"),
    :filters => PropDict(
        :sg => PropDict(
            :min_tot_intersect => 100u"ns", :max_tot_intersect => 1000u"ns",
            :min_threshold => -50.0, :max_threshold => 50.0, :n_σ_threshold => 5.0,
            :min_dc_threshold => -50.0, :max_dc_threshold => 50.0, :n_σ_dc_threshold => 5.0,
        ),
        :trap => PropDict(
            :rt => 200u"ns", :ft => 100u"ns", :pz_tau => 1000u"µs",
            :min_tot_intersect => 100u"ns", :max_tot_intersect => 1000u"ns",
            :min_threshold => -50.0, :max_threshold => 50.0, :n_σ_threshold => 5.0,
            :min_dc_threshold => -50.0, :max_dc_threshold => 50.0, :n_σ_dc_threshold => 5.0,
        ),
    ),
)
pars_sipm = PropDict(:sg => PropDict(:wl => 0.3u"µs"))

sipm_out = dsp_sipm(data, sipm_cfg, pars_sipm)
(length.(sipm_out.trig_pos_sg), round.(u"µs", sipm_out.trig_pos_sg[2], digits = 2))
```

`dsp_sipm_compressed` returns the same table when it is given the same data:

```@example sipm_pmt
sipm_out_bit_drop = dsp_sipm_compressed(data, sipm_cfg, pars_sipm)
sipm_out_bit_drop.trig_pos_sg == sipm_out.trig_pos_sg
```

The following figure shows the second waveform with the positions of the SG triggers.

```@example sipm_pmt
using CairoMakie
using RadiationDetectorDSP

wvf_deriv = SavitzkyGolayFilter(0.3u"µs", 2, 1)(wvfs[2])

fig = Figure(size = (700, 400))
ax1 = Axis(fig[1, 1], ylabel = "Signal (ADC counts)")
lines!(ax1, ustrip.(u"µs", wvfs[2].time), wvfs[2].signal)
vlines!(ax1, ustrip.(u"µs", sipm_out.trig_pos_sg[2]), color = :red, linestyle = :dash)
ax2 = Axis(fig[2, 1], xlabel = "Time (µs)", ylabel = "SG derivative")
lines!(ax2, ustrip.(u"µs", wvf_deriv.time), wvf_deriv.signal)
hlines!(ax2, sipm_out.threshold_sg[2], color = :gray, linestyle = :dash)
vlines!(ax2, ustrip.(u"µs", sipm_out.trig_pos_sg[2]), color = :red, linestyle = :dash)
linkxaxes!(ax1, ax2)
fig
```

## SiPM Parameter Scans

The window length of the Savitzky-Golay filter and the trigger thresholds are determined from data by two scan functions.
They take the waveforms of a SiPM channel with a baseline of about zero.
The sequence is:

1. [`dsp_sg_sipm_optimization_compressed`](@ref) scans a grid of window lengths and returns, for each window length, the maxima of all triggers found in the derivative and the threshold that was used.
2. A fit of the single-photoelectron peak in the distribution of maxima as a function of the window length selects the window length. This fit is not part of LegendDSP.jl; it is done with `fit_sipm_wl` of LegendSpecFits.jl.
3. [`dsp_sg_sipm_thresholds_compressed`](@ref) returns the samples of the derivative and of the integrated waveform at that window length, so that the distribution of noise values can be fitted to set the thresholds of the chain.

The result of the scan is stored as `pars_optimization` for [`dsp_sipm`](@ref).

[`dsp_sg_sipm_optimization_compressed`](@ref)`(wvfs, dsp_config, optimization_config)` reads from `dsp_config` the keys `min_tot_intersect`, `max_tot_intersect`, `n_σ_threshold`, and `sg_flt_degree`, which are at the top level here and not below `filters.sg`.
From `optimization_config` it reads `e_grid_wl`, a collection of window lengths (times), `threshold.min_cut` and `threshold.max_cut`, the bounds of the samples that enter the noise estimate, and `threshold.n_wvfs`, the number of waveforms (counted from the first) used for the noise estimate.
For each window length, it

1. filters all waveforms with the Savitzky-Golay filter of that window length (polynomial degree `sg_flt_degree`, derivative order 1),
2. computes the standard deviation ([`thresholdstats`](@ref)) of the samples of the first `threshold.n_wvfs` filtered waveforms that lie between `min_cut` and `max_cut` and multiplies it by `n_σ_threshold`, and
3. finds the maxima of the triggers above this threshold in all waveforms with [`IntersectMaximum`](@ref).

It returns a `Table` with the columns `trig_max_grid`, a `VectorOfVectors` with the trigger maxima of all waveforms for each window length, and `thresholds_grid`, the threshold for each window length.

The method `dsp_sg_sipm_optimization_compressed(n_max_wvfs, wvfs, dsp_config, optimization_config)` is intended to process the waveforms in chunks of at most `n_max_wvfs` waveforms and to combine the chunks, with the smallest threshold of the chunks for each window length.
It throws a `MethodError` when it is called.

[`dsp_sg_sipm_thresholds_compressed`](@ref)`(wvfs, sg_window_length, config)` reads `config.sg_flt_degree`.
It filters the waveforms with the derivative Savitzky-Golay filter, integrates the result, and returns a `Table` with the columns `bsl_deriv` (the samples of all filtered waveforms in one vector), `bsl` (the samples of the integrated waveforms in one vector), and `bsl_flipped` (`-bsl`).
All waveforms must have the same number of samples.

```@example sipm_pmt
scan_cfg = PropDict(
    :sg_flt_degree => 2, :min_tot_intersect => 100u"ns", :max_tot_intersect => 1000u"ns", :n_σ_threshold => 5.0,
)
scan_opt = PropDict(
    :e_grid_wl => range(0.1u"µs", 0.5u"µs", step = 0.1u"µs"),
    :threshold => PropDict(:min_cut => -50.0, :max_cut => 50.0, :n_wvfs => 3),
)
wvfs_zero = shift_waveform.(wvfs, -20.0)

grid = dsp_sg_sipm_optimization_compressed(wvfs_zero, scan_cfg, scan_opt)
(length.(grid.trig_max_grid), round.(grid.thresholds_grid, digits = 3))
```

The threshold falls with the window length because longer windows average more noise.

```@example sipm_pmt
thresholds = dsp_sg_sipm_thresholds_compressed(wvfs_zero, 0.3u"µs", PropDict(:sg_flt_degree => 2))
(length(thresholds.bsl_deriv), propertynames(thresholds))
```

## PMT Chain

[`dsp_pmts`](@ref)`(data, config)` processes the waveforms of one PMT channel.
The input table needs the columns `waveform`, `timestamp`, `channel`, `eventnumber`, and `daqenergy`.
There is no separate function for compressed data.

The steps are:

1. [`TimeAxisFilter`](@ref)`(time_axis_step_length)` sets the sampling period of the waveforms.
2. The baseline mean, standard deviation, and slope are determined between `baseline_window_start` and `baseline_window_end` with `signalstats`, and the mean is subtracted.
3. The maximum and minimum of the baseline-subtracted waveforms and their times are determined ([`extremestats`](@ref)); these are the `raw_` columns.
4. [`IntersectMaximum`](@ref)`(min_tot_intersect, max_tot_intersect)` finds all crossings of the fixed threshold `intersect_threshold` and the maximum after each crossing.
5. [`saturation`](@ref) counts the samples at `saturation_limit_low` and `saturation_limit_high` in the waveform before the baseline subtraction.
6. The baseline-subtracted waveforms are smoothed. If `wsg_weight` is 0, a `SavitzkyGolayFilter(wsg_window_length, wsg_flt_degree, 0)` of RadiationDetectorDSP.jl is used, otherwise a [`WeightedSavitzkyGolayFilter`](@ref) with `weightType = wsg_weight`. The maximum and minimum of the smoothed waveforms and their times are the columns without the `raw_` prefix.

The configuration keys are:

| Key | Meaning |
|:---|:---|
| `time_axis_step_length` | sampling period assigned to the waveforms (time) |
| `baseline_window_start`, `baseline_window_end` | limits of the baseline window (times) |
| `min_tot_intersect`, `max_tot_intersect` | minimum time over threshold and window for the maximum (times of the same unit) |
| `intersect_threshold` | trigger threshold in units of the baseline-subtracted signal |
| `wsg_window_length`, `wsg_flt_degree` | window length (time) and polynomial degree of the smoothing filter |
| `wsg_weight` | weight function of the smoothing filter; 0 selects the unweighted filter, 1 to 4 select the weights of [`WeightedSavitzkyGolayFilter`](@ref) |
| `saturation_limit_high`, `saturation_limit_low` | digitizer limits for the saturation count |

The output table has the following columns.

| Column | Content |
|:---|:---|
| `timestamp`, `eventID_fadc`, `e_fc`, `channel` | the input columns `timestamp`, `eventnumber`, `daqenergy`, and `channel` |
| `raw_pulse_height`, `raw_pulse_low` | maximum and minimum of the baseline-subtracted waveform |
| `raw_t0_hi`, `raw_t0_low` | times of the maximum and minimum |
| `trig_max`, `trig_t`, `trig_mult` | maxima and times of the triggers (vectors per event) and their number |
| `sat_low`, `sat_high` | number of samples at the lower and upper digitizer limit |
| `pulse_height`, `pulse_low`, `t0_hi`, `t0_low` | maximum, minimum and their times after smoothing |
| `bl_mean`, `bl_sigma`, `bl_slope` | baseline statistics |

```@example sipm_pmt
pmt_cfg = PropDict(
    :time_axis_step_length => 4u"ns",
    :baseline_window_start => 0u"µs", :baseline_window_end => 4u"µs",
    :min_tot_intersect => 40u"ns", :max_tot_intersect => 400u"ns",
    :intersect_threshold => 3.0,
    :wsg_window_length => 100u"ns", :wsg_flt_degree => 2, :wsg_weight => 2,
    :saturation_limit_high => 16383, :saturation_limit_low => 0,
)

t_pmt = range(0.0u"µs", step = 0.004u"µs", length = 2000)
wvfs_pmt = ArrayOfRDWaveforms([sipm_waveform(rng, t_pmt; pulses = p, tau = 0.2, noise = 0.3) for p in pulses])
data_pmt = Table(
    waveform = wvfs_pmt, timestamp = [1.0, 2.0, 3.0], eventnumber = [1, 2, 3],
    daqenergy = [100, 200, 300], channel = [1, 1, 1],
)

pmt_out = dsp_pmts(data_pmt, pmt_cfg)
(pmt_out.trig_mult, round.(pmt_out.pulse_height, digits = 2))
```

The second pulse of the second waveform starts at 15 µs, after the end of the 8 µs long PMT waveform, so only one trigger is found.
