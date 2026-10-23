# HPGe DSP Chains

A DSP chain takes the raw waveform table of one detector channel and returns a table with one row per event and one column per extracted parameter. The HPGe chain [`dsp_icpc`](@ref) and its variant for compressed raw data [`dsp_icpc_compressed`](@ref) compute the baseline, the pole-zero (PZ) corrected tail, timing points, energy estimators from several filters, the maximum current (A), the charge-drift parameter (Q-drift), late-charge parameters (LQ), pile-up and saturation flags, and an optional machine-learning quality-cut label. The pulser and auxiliary channels use the shorter chains [`dsp_puls`](@ref), [`dsp_puls_compressed`](@ref), [`dsp_aux`](@ref) and [`dsp_aux_compressed`](@ref).

The layout of the input tables and the [`DSPConfig`](@ref) are described in [Waveforms and Configuration](@ref). The individual filters are described in [Filters](@ref), and the functions that compute single parameters from waveforms in [Waveform Parameters](@ref).

## Calling a Chain

```julia
out = dsp_icpc(data, config, τ, pars_filter; f_evaluate_qc = missing)
out = dsp_icpc_compressed(data, config, τ, pars_filter; f_evaluate_qc = missing)
```

| Argument | Meaning |
|:---------|:--------|
| `data` | A `TypedTables.Table` with the raw columns listed below |
| `config` | A [`DSPConfig`](@ref) |
| `τ` | The decay time of the preamplifier signal as a Unitful time (for example `450.0u"µs"`). It sets the pole-zero correction. How `τ` is determined is described in [Decay Time](@ref) |
| `pars_filter` | A `PropDict` with the optimized filter parameters of the channel (`trap`, `cusp`, `zac` with `rt` and `ft`; `sg` and `mwa` with `wl`). Missing entries are taken from `config.default_flt_param` |
| `f_evaluate_qc` | Keyword. Either `missing` (default) or a function that evaluates a trained classifier, see [Quality-Cut Label](@ref) |

The input columns are

* `dsp_icpc`: `waveform`, `baseline`, `timestamp`, `eventnumber`, `daqenergy`.
* `dsp_icpc_compressed`: `waveform_presummed`, `waveform_windowed`, `presum_rate`, `baseline`, `timestamp`, `eventnumber`, `daqenergy`, `t_sat_lo`, `t_sat_hi`, `deadtime`.

The columns `baseline`, `timestamp`, `eventnumber` and `daqenergy` are the values reported by the FADC. They are copied to the output as `blfc`, `timestamp`, `eventID_fadc` and `e_fc`. Encoded columns are decoded first.

## Processing Steps

The two variants run the same sequence. In the compressed variant the steps that need the full trace use `waveform_presummed` (the "presummed waveform") and the steps that need fine time resolution use `waveform_windowed` (the "windowed waveform"); see [Uncompressed and Compressed Input](@ref) in [Waveforms and Configuration](@ref).

1. **Saturation.** The samples equal to the lowest ADC value (0) and to the highest value (`2^fc_bit_depth - fc_bit_depth`, multiplied by `presum_rate` for presummed waveforms) are counted, together with the longest run of consecutive saturated samples at each end ([`saturation`](@ref)).
2. **Baseline.** A linear fit in `bl_window` gives mean, standard deviation, slope and offset of the baseline (`signalstats`). The mean is subtracted from the waveforms.
3. **Quality-cut label.** If `f_evaluate_qc` is given, it is evaluated on the baseline-subtracted waveform (see below). Otherwise the label is `-1`.
4. **Extrema.** The maximum and minimum of the baseline-subtracted waveform, before the PZ correction.
5. **Tail.** [`tailstats`](@ref) fits a straight line to the logarithm of the baseline-subtracted waveform in `tail_window`. The result is the tail decay time (`tail_τ`) and the mean and standard deviation of the logarithm. See [Decay Time](@ref).
6. **Pole-zero correction.** The waveform is deconvolved with `InvCRFilter(τ)`, which turns the exponential tail into a flat top. Steps 7 to 17 work on the PZ-corrected waveform. The tail window is evaluated again (`signalstats`) to check the correction.
7. **t0.** [`get_t0`](@ref) applies a fast asymmetric trapezoid with the lengths `kwargs_pars.t0_flt_pars` and finds the first time the output stays above `t0_threshold` for at least `kwargs_pars.t0_mintot`.
8. **Fixed-length energy filters.** Trapezoidal filters with (rise, flat-top) of 10 µs and 4 µs (`e_10410`), 5 µs and 3 µs (`e_535`), 3 µs and 1 µs (`e_313`) and 10 µs and 10 µs (`e_101010`). The energy is the maximum of the filter output.
9. **Rise-time points.** `t10`, `t50`, `t80`, `t90`, `t99` are the times at which the PZ-corrected waveform first stays for `kwargs_pars.tx_mintot` above 10 %, 50 %, 80 %, 90 % and 99 % of the maximum of the 10 µs/10 µs trapezoid ([`get_threshold`](@ref)). The drift time is `t90 - t0`.
10. **Maximum of the simple derivative.** [`DerivativeFilter`](@ref) gives the current without smoothing. The time of its maximum is `t_amax`.
11. **Q-drift.** [`get_qdrift`](@ref) integrates the waveform and compares two consecutive integration windows that start at `t0`, delimited by the first and last element of `qdrift_int_length`.
12. **LQ.** For every `[early, gap, late]` triple of `lq_windows` and for each of the start points `t80`, `t90` and `t_amax`, [`get_lq`](@ref) subtracts the charge in the early window from the charge in the late window.
13. **Optimized energy filters.** Trapezoidal, CUSP and ZAC filters with the optimized (or default) parameters are applied. The energy is not the filter maximum but the filter output at a fixed pick-off time after `t50`, read with a polynomial signal estimator (`kwargs_pars.sig_interpolation_order` and `sig_interpolation_length`). The pick-off time is `t50 + rt + ft/2` for the trapezoid and `t50 + flt_length/2` for CUSP and ZAC. The maximum of each filter output and its time are stored as well.
14. **Current (A).** The maximum in `current_window` of three versions of the current: the Savitzky-Golay derivative with window `sg.wl` and degree `sg_flt_degree` (`a_sg`), the simple derivative (`a_raw`), and the simple derivative smoothed by `a_mwa_passes` alternating passes of a moving average of window `mwa.wl` (`a_mwa`). Additional columns are written for the fixed windows `a_sg_windows` and `a_mwa_windows`. The maxima are interpolated quadratically around the maximum sample.
15. **Pile-up.** [`get_triggers`](@ref) counts the intersections of a smoothed signal with a threshold of `inTraceCut_std_threshold` times its standard deviation in `bl_window`, with a minimum time over threshold of `kwargs_pars.intrace_mintot`. `dsp_icpc` applies it to the Savitzky-Golay derivative; `dsp_icpc_compressed` applies it to the presummed waveform after a 2 µs/1 µs trapezoidal filter with a minimum time over threshold of 2 µs.
16. **Current rise.** `t50_current` is the time at which the SG derivative of the whole waveform first stays above half of its maximum.
17. **Inverted waveform.** The PZ-corrected waveform is multiplied by `-1`. The maxima of the 10 µs/4 µs and 3 µs/1 µs trapezoids (`e_10410_inv`, `e_313_inv`) and the t0 of the inverted waveform (`t0_inv`) are stored; they flag events with a negative-going component after the pulse.

In the compressed variant the pole-zero corrected presummed waveform is additionally evaluated in the two windows `auxpz1_window` and `auxpz2_window`.

## A Synthetic Example

The following example runs both variants on synthetic waveforms: a linear rise of 1 µs followed by an exponential tail, white noise, and baseline 1000. The configuration is the generic example from [Waveforms and Configuration](@ref).

```@example dsp_chains
using LegendDSP
using PropDicts, Unitful, StableRNGs, TypedTables, RadiationDetectorSignals, ArraysOfArrays

grid(start, stop, step) = PropDict(:start => start, :stop => stop, :step => step)
rt_ft_grid = PropDict(:rt => grid(1.0u"µs", 16.0u"µs", 0.5u"µs"), :ft => grid(1.0u"µs", 4.0u"µs", 0.2u"µs"))

config = DSPConfig(PropDict(
    :enc_pickoff_trap => 40.0u"µs", :enc_pickoff_zac => 41.0u"µs", :enc_pickoff_cusp => 41.0u"µs",
    :flt_length_cusp => 38.0u"µs", :flt_length_zac => 38.0u"µs",
    :bl_window => PropDict(:min => 0.0u"µs", :max => 39.0u"µs"),
    :tail_window => PropDict(:min => 70.0u"µs", :max => 110.0u"µs"),
    :current_window => PropDict(:min => 43.0u"µs", :max => 62.0u"µs"),
    :auxpz1_window => PropDict(:min => 70.0u"µs", :max => 90.0u"µs"),
    :auxpz2_window => PropDict(:min => 90.0u"µs", :max => 110.0u"µs"),
    :t0_threshold => 4.0, :inTraceCut_std_threshold => 5, :sg_flt_degree => 3,
    :qdrift_int_length => [2.5u"µs", 5.0u"µs"],
    :lq_windows => [[2.5, 0.0, 2.5]u"µs"],
    :e_grid_trap => rt_ft_grid, :e_grid_zac => rt_ft_grid, :e_grid_cusp => rt_ft_grid,
    :a_grid_wl_sg => grid(30.0u"ns", 350.0u"ns", 32.0u"ns"),
    :a_mwa_windows => [48.0u"ns", 288.0u"ns"], :a_mwa_passes => 3,
    :a_sg_windows => [60.0u"ns", 100.0u"ns"],
    :flt_defaults => PropDict(
        :sg => 100.0u"ns", :mwa => 96.0u"ns",
        :trap => PropDict(:rt => 5.0u"µs", :ft => 2.5u"µs"),
        :zac => PropDict(:rt => 5.0u"µs", :ft => 2.5u"µs"),
        :cusp => PropDict(:rt => 5.0u"µs", :ft => 2.5u"µs"),
    ),
    :kwargs_pars => PropDict(
        :fc_bit_depth => 16,
        :t0_flt_pars => [40.0u"ns", 100.0u"ns", 2000.0u"ns"],
        :t0_mintot => 1500.0u"ns", :tx_mintot => 32.0u"ns", :intrace_mintot => 100.0u"ns",
        :int_interpolation_order => 3, :int_interpolation_length => 100.0u"ns",
        :sig_interpolation_order => 3, :sig_interpolation_length => 700.0u"ns",
    ),
))
nothing # hide
```

The waveform generator and the raw table of the uncompressed layout:

```@example dsp_chains
rng = StableRNG(1)
τ_true = 450.0u"µs"

function synthetic_waveform(amplitude; n = 8192, dt = 16.0u"ns", t_rise = 50.0u"µs", rise = 1.0u"µs")
    t = range(0.0u"ns", step = dt, length = n)
    x = ustrip.(u"µs", t) .- ustrip(u"µs", t_rise)
    shape = map(x) do xi
        xi < 0 ? 0.0 : xi < ustrip(u"µs", rise) ? xi / ustrip(u"µs", rise) :
            exp(-(xi - ustrip(u"µs", rise)) / ustrip(u"µs", τ_true))
    end
    RDWaveform(t, 1000 .+ amplitude .* shape .+ 3 .* randn(rng, n))
end

amplitudes = [4500.0, 5000.0, 5500.0, 6000.0]
N = length(amplitudes)
wvfs = ArrayOfRDWaveforms([synthetic_waveform(a) for a in amplitudes])

raw = Table(
    waveform = wvfs, baseline = fill(1000.0f0, N), timestamp = zeros(UInt64, N),
    eventnumber = UInt32.(1:N), daqenergy = zeros(UInt16, N),
)
nothing # hide
```

Running the chain with an empty `pars_filter` uses the default filter parameters of the configuration:

```@example dsp_chains
out = dsp_icpc(raw, config, τ_true, PropDict())
length(columnnames(out))
```

A selection of the output columns. The three energy estimators are close to the amplitude of the baseline-subtracted signal:

```@example dsp_chains
Table(amplitude = amplitudes, e_trap = out.e_trap, e_cusp = out.e_cusp, e_zac = out.e_zac)
```

The timing points are in µs, and the drift time is `t90 - t0`. The tail decay time before the pole-zero correction recovers the decay time that was used to generate the waveforms:

```@example dsp_chains
Table(t0 = out.t0, t50 = out.t50, t90 = out.t90, drift_time = out.drift_time,
      tail_τ = uconvert.(u"µs", out.tail_τ))
```

Optimized filter parameters override the defaults entry by entry:

```@example dsp_chains
pars_filter = PropDict(:trap => PropDict(:rt => 8.0u"µs", :ft => 2.0u"µs"))
out_opt = dsp_icpc(raw, config, τ_true, pars_filter)
(default = out.e_trap, optimized = out_opt.e_trap)
```

The compressed variant needs a presummed and a windowed waveform per event. The following lines build them from the full-rate waveforms by summing groups of four samples and by cutting a window of about 30 µs around the rise:

```@example dsp_chains
presum_rate = 4

function presummed(w, r)
    s = reshape(w.signal, r, :)
    RDWaveform(range(first(w.time), step = r * step(w.time), length = size(s, 2)), vec(sum(s, dims = 1)))
end
windowed(w, idx) = RDWaveform(w.time[idx], w.signal[idx])

raw_compressed = Table(
    waveform_presummed = ArrayOfRDWaveforms([presummed(w, presum_rate) for w in wvfs]),
    waveform_windowed = ArrayOfRDWaveforms([windowed(w, 2500:4400) for w in wvfs]),
    presum_rate = fill(UInt16(presum_rate), N),
    baseline = fill(1000.0f0, N), timestamp = zeros(UInt64, N),
    eventnumber = UInt32.(1:N), daqenergy = zeros(UInt16, N),
    t_sat_lo = zeros(UInt16, N), t_sat_hi = zeros(UInt16, N), deadtime = zeros(UInt16, N),
)

out_c = dsp_icpc_compressed(raw_compressed, config, τ_true, PropDict())
Table(amplitude = amplitudes, e_trap = out_c.e_trap, e_trap_over_rate = out_c.e_trap ./ presum_rate,
      t50 = out_c.t50, t50_pre = out_c.t50_pre, trigger_multiplicity = out_c.trigger_multiplicity)
```

The energies of the compressed variant are in units of presummed ADC counts, `presum_rate` times larger than the ADC amplitude. The time points `t0` to `t99` come from the windowed waveform and `t50_pre` from the presummed waveform; the latter has a coarser time resolution.

Pulser and auxiliary channels are processed by the shorter chains:

```@example dsp_chains
puls = dsp_puls(raw, config)
aux = dsp_aux(raw, config)
(puls_columns = columnnames(puls), aux_columns = columnnames(aux))
```

## Output Columns

Times are Unitful quantities; their unit is given below. Columns without a unit are in ADC counts (for the compressed variant see the note on presummed values above), counts, or dimensionless. Entries marked "c" exist only in `dsp_icpc_compressed`; entries marked "u" exist only in `dsp_icpc`.

### FADC Values

| Column | Meaning |
|:-------|:--------|
| `blfc` | Baseline reported by the FADC (input `baseline`) |
| `timestamp` | Event timestamp (input) |
| `eventID_fadc` | Event number (input `eventnumber`) |
| `e_fc` | Energy reported by the FADC (input `daqenergy`) |
| `deadtime` (c) | Input `deadtime`, copied |
| `t_sat_lo`, `t_sat_hi` (c) | Input columns, copied |

### Baseline and Tail

| Column | Meaning |
|:-------|:--------|
| `blmean`, `blsigma` | Mean and standard deviation of the samples in `bl_window` |
| `blslope`, `bloffset` | Slope (1/time unit of the time axis) and offset of the linear fit in `bl_window` |
| `bl_slope_sigma` (c) | Standard deviation of the residual after the linear fit in `bl_window` |
| `tail_τ` | Decay time from the exponential fit of the tail, before the PZ correction (time unit of the time axis) |
| `tail_mean`, `tail_sigma` | Mean and standard deviation of the natural logarithm of the samples in `tail_window`, before the PZ correction |
| `tailmean`, `tailsigma`, `tailslope`, `tailoffset` | Mean, standard deviation, slope and offset of the PZ-corrected waveform in `tail_window` |
| `auxpz1_mean`, `auxpz1_sigma`, `auxpz1_slope_sigma`, `auxpz2_mean`, `auxpz2_sigma`, `auxpz2_slope_sigma` (c) | Mean, standard deviation and residual standard deviation of the PZ-corrected presummed waveform in `auxpz1_window` and `auxpz2_window` |

A flat PZ-corrected tail has `tailslope` close to zero. A negative slope indicates that `τ` is too large, a positive slope that it is too small (see [Decay Time](@ref)).

### Extrema

| Column | Meaning |
|:-------|:--------|
| `e_max`, `e_min` | Maximum and minimum of the baseline-subtracted waveform before the PZ correction. In the compressed variant `e_max` is taken from the presummed waveform and `e_min` from the windowed waveform |
| `t_max`, `t_min` (c) | Times of `e_max` and `e_min` |
| `e_max_pre`, `e_min_pre`, `t_max_pre`, `t_min_pre` (c) | Maxima, minima and their times of the presummed waveform |

### Timing

| Column | Meaning |
|:-------|:--------|
| `t0` | Start of the signal rise from the fast asymmetric trapezoid, in µs; `0` if no crossing is found |
| `t10`, `t50`, `t80`, `t90`, `t99` | Times of 10 %, 50 %, 80 %, 90 %, 99 % of the 10 µs/10 µs trapezoid maximum, in µs; `0` if no crossing is found |
| `t50_pre` (c) | Same as `t50`, found on the presummed waveform |
| `drift_time` | `t90 - t0` in ns |
| `t50_current` | Time at which the SG-smoothed current first reaches half of its maximum, in µs |
| `t_amax` | Time of the maximum of the simple derivative, quadratically interpolated, in the unit of the time axis |
| `t_a_sg`, `t_a_mwa` | Times of the maxima `a_sg` and `a_mwa`, in the unit of the time axis |
| `t0_inv` | t0 of the inverted PZ-corrected waveform, in µs |

### Energy

| Column | Meaning |
|:-------|:--------|
| `e_10410`, `e_535`, `e_313`, `e_101010` | Maximum of the trapezoid (rise, flat-top) of 10/4, 5/3, 3/1 and 10/10 µs |
| `t_10410`, `t_535`, `t_313`, `t_101010` (c) | Times of those maxima |
| `e_trap`, `e_cusp`, `e_zac` | Energy from the trapezoidal, CUSP and ZAC filter with the optimized or default parameters, read at the pick-off time |
| `e_trap_max`, `e_cusp_max`, `e_zac_max` | Maximum of the output of those filters |
| `t_trap_max`, `t_cusp_max`, `t_zac_max` | Times of those maxima, in the unit of the time axis |
| `e_10410_inv`, `e_313_inv` | Maximum of the 10/4 and 3/1 µs trapezoid on the inverted PZ-corrected waveform |

For the CUSP and ZAC filters `rt` is the shaping time σ and `ft` the length of the flat top; the total length of the filter is `flt_length_cusp` or `flt_length_zac`. The decay-time constant of these filters is set to a very long value so that the filters do not apply a CR correction of their own.

### Current, Q-drift and LQ

| Column | Meaning |
|:-------|:--------|
| `a_sg` | Maximum of the Savitzky-Golay derivative in `current_window`, in ADC counts per sample |
| `a_raw` | Maximum of the simple derivative in `current_window`, in ADC counts per sample |
| `a_mwa` | Maximum of the moving-window-averaged derivative in `current_window`, in ADC counts per sample |
| `a_sg_<ns>`, `a_mwa_<ns>` | The same for the fixed window lengths in `a_sg_windows` and `a_mwa_windows`, with the window length in ns in the name (`a_sg_60`) |
| `qdrift` | Difference of the integrated charge of two consecutive windows starting at `t0` |
| `lq80_<e><g><l>`, `lq90_<e><g><l>`, `lqamax_<e><g><l>` | LQ with start point `t80`, `t90` or `t_amax`. The name encodes the early, gap and late window lengths in units of 0.1 µs with two digits each: `lq80_250025` is the LQ with 2.5 µs early window, no gap and 2.5 µs late window. One column per start point and per triple in `lq_windows` |

The ratio of `a_sg` (or `a_mwa`) to an energy estimator is the A/E parameter used for pulse-shape discrimination; see the [A/E Pulse-Shape Parameter](@ref) tutorial.

### Pile-up and Saturation

| Column | Meaning |
|:-------|:--------|
| `inTrace_intersect`, `inTrace_n` (u) | Position (unit of the time axis) and number of triggers found in the SG-smoothed current when the threshold is `inTraceCut_std_threshold` times the baseline standard deviation of that signal |
| `trigger_position`, `trigger_multiplicity` (c) | Position and number of triggers found in a 2 µs/1 µs trapezoid of the presummed waveform with the same threshold definition and a minimum time over threshold of 2 µs |
| `n_sat_low`, `n_sat_high` | Number of samples at the lowest and highest ADC value |
| `n_sat_low_cons`, `n_sat_high_cons` | Longest run of consecutive samples at those values |
| `qc_label` | Label of the machine-learning classifier, `-1` if none is used |

For the presummed waveforms of the compressed variant the saturation counts are made on the presummed waveform, with the highest value scaled by `presum_rate`.

## Quality-Cut Label

`f_evaluate_qc` is `missing` by default and the chain then sets `qc_label = -1` for every event without running any classifier. A function can be passed instead. The chain prepares its input in the following way: the baseline-subtracted waveform (before the PZ correction) is smoothed by repeated [`HaarAveragingFilter`](@ref) passes (5 passes for `dsp_icpc`, 2 passes for `dsp_icpc_compressed`; see [`get_qc_classifier`](@ref) and [`get_qc_classifier_compressed`](@ref)), normalized to a maximum absolute value of one, and collected in a matrix with one column per event. The function receives this matrix and returns a tuple `(labels, decision_values)`. The labels are stored in `qc_label` as integers. [Machine-Learning Quality Cuts](@ref) describes how to train such a function with [`get_qc_ml_func`](@ref).

```@example dsp_chains
# A stand-in for a trained classifier: label every event with 1.
toy_classifier(x) = (ones(size(x, 2)), zeros(size(x, 2)))

out_qc = dsp_icpc(raw, config, τ_true, PropDict(); f_evaluate_qc = toy_classifier)
out_qc.qc_label |> collect
```

## Pulser and Auxiliary Chains

| Function | Reads | Output columns |
|:---------|:------|:---------------|
| [`dsp_puls`](@ref) | `waveform` | `blmean`, `blsigma`, `blslope`, `bloffset` (baseline in `bl_window`), `t50` (µs, time of half of the baseline-subtracted maximum), `e_max` (maximum), `e_10410` (maximum of the 10 µs/4 µs trapezoid), and the FADC values |
| [`dsp_puls_compressed`](@ref) | `waveform_presummed` | the same columns, computed on the presummed waveform |
| [`dsp_aux`](@ref) | `waveform` | `e_max` (maximum minus the FADC baseline `baseline`), `t_max` (time of the maximum), and the FADC values |
| [`dsp_aux_compressed`](@ref) | `waveform_presummed` | `e_max` (maximum minus `baseline * presum_rate`), `t_max`, and the FADC values |

The pulser chains read `bl_window` of the configuration; the auxiliary chains do not read the configuration. `dsp_aux` and `dsp_aux_compressed` are meant for channels whose waveforms are spikes, such as veto channels.
