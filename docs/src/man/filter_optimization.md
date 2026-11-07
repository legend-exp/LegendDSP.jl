# Filter Optimization

The energy filters (trapezoidal, CUSP, ZAC) and the current filters (Savitzky-Golay, moving window average) of the HPGe chains have free parameters: the rise time `rt` and flat-top time `ft` of an energy filter and the window length `wl` of a current filter. The optimal values depend on the detector and its electronics. LegendDSP.jl provides functions that evaluate a filter on a grid of parameter values for a set of waveforms of one detector. They return the filter output for every grid point and every waveform. The statistical evaluation of these grids (noise width, energy resolution, A/E separation) is done by the fit functions of LegendSpecFits.jl; see the [LegendSpecFits.jl documentation](https://legend-exp.github.io/LegendSpecFits.jl/stable/). The optimal values enter the DSP chains as the argument `pars_filter` ([HPGe DSP Chains](@ref)).

The energy filter scans share one preparation step, [`dsp_flt_optimization_prepare`](@ref): encoded waveforms are decoded, the baseline mean in `config.bl_window` is subtracted, the waveforms are corrected for the pole-zero with `InvCRFilter(τ)` (see [Decay Time](@ref)), and `t50` is determined. The result is a [`FltOptPrep`](@ref). The scans therefore need the decay time `τ` of the detector. The grids and the other settings come from the [`DSPConfig`](@ref). The scans that take waveforms and `τ` prepare them internally on every call; the [prepared workflow](@ref "Prepared Workflow") does the preparation once.

## Overview

| Function | Scanned parameter | Grid in `config` | Result |
|:---------|:------------------|:-----------------|:-------|
| [`dsp_trap_rt_optimization`](@ref), [`dsp_cusp_rt_optimization`](@ref), [`dsp_zac_rt_optimization`](@ref) | rise time `rt` at fixed flat-top time (keyword `ft`, default 2 µs); take a [`FltOptPrep`](@ref) or waveforms and `τ` | `e_grid_rt_trap`, `e_grid_rt_cusp`, `e_grid_rt_zac` | `Matrix{Float64}`, rows: grid points, columns: waveforms |
| [`dsp_trap_ft_optimization`](@ref), [`dsp_cusp_ft_optimization`](@ref), [`dsp_zac_ft_optimization`](@ref) | flat-top time `ft` at fixed rise time `rt` (argument); take waveforms and `τ` | `e_grid_ft_trap`, `e_grid_ft_cusp`, `e_grid_ft_zac` | `Matrix{Union{Missing,Float32}}`, rows: grid points, columns: waveforms (every entry is set; no `missing` is stored) |
| [`dsp_flt_optimization_prepare`](@ref), [`dsp_energy_rtft`](@ref), [`realized_flattop`](@ref) | single points of the (`rt`, `ft`) plane; see [Prepared Workflow](@ref) | none | `FltOptPrep`, energies of all waveforms at one point, flat-top time that a filter applies |
| [`dsp_qc_flt_optimization`](@ref), [`dsp_qc_flt_optimization_compressed`](@ref) | none | none | `Table` with the energy at the default trapezoid, baseline mean and slope, `t50`, and `qc_label` |
| [`dsp_dplms_optimization`](@ref) | weights and flat-top widening of the DPLMS filter, see [DPLMS Filter](@ref); takes waveforms, `τ` and a drift-time proxy | none (`grid` argument, [`dplms_default_grid`](@ref)) | `NamedTuple` with the energies of every grid point, the filters, and the matrices |
| [`dsp_qdrift_flt_optimization`](@ref) | none | `qdrift_int_length` | Vector of Q-drift values |
| [`dsp_sg_optimization`](@ref), [`dsp_sg_optimization_compressed`](@ref) | window length of the Savitzky-Golay current filter | `a_grid_wl_sg` | `Table` with A/E for every window length |
| [`dsp_mwa_optimization`](@ref), [`dsp_mwa_optimization_compressed`](@ref) | window length of the moving-window-average current filter | `a_grid_wl_mwa` | `Table` with A/E for every window length |

There are no `_compressed` variants of the rise-time and flat-top scans. They work on a single waveform per event, and the dataflow passes them the presummed waveforms. In the rise-time scans with waveforms and `τ`, the keyword `ft` and `τ` must have the same number type (for example both `Float64`); the same holds for `rt` and `τ` in the flat-top scans.

## Energy Filter Scans

### Rise Time

`dsp_{trap,cusp,zac}_rt_optimization(wvfs, config, τ; ft)` applies the filter with every rise time in the grid to all waveforms and reads the filter output at the fixed time `config.enc_pickoff_{trap,cusp,zac}`. This pick-off time lies in the part of the trace before the signal rise, so the values are samples of the filter-output noise and no signal. Their spread over the waveforms, as a function of the rise time, is the electronic noise (ENC) of the filter. The function returns the values; the width is determined by the fit. For the CUSP and ZAC filters the grid value is the shaping time σ, and `ft` is the length of the flat top; the total filter length is `config.flt_length_cusp` or `config.flt_length_zac`.

### Flat-Top Time

`dsp_{trap,cusp,zac}_ft_optimization(wvfs, config, τ, rt)` applies the filter with the given rise time `rt` and every flat-top time in the grid. The energy is read at a pick-off time that depends on the filter: `t50 + rt + ft/2` for the trapezoid and `t50 + flt_length/2` for CUSP and ZAC, where `t50` is the time at which the PZ-corrected waveform reaches half of its maximum. The result is the energy of every waveform for every flat-top time, stored as `Float32`. An energy is `NaN` if the pick-off time lies outside the time axis of the filtered waveform. Row `i` equals `dsp_energy_rtft(prep, config, rt, ft_i, flt)` for the prepared waveforms, up to the conversion to `Float32`. A flat-top time shorter than the rise time of the signal gives a lower energy (ballistic deficit). For a calibration line with many events, the width of the line as a function of the flat-top time shows the flat-top time that gives the best resolution. The Q-drift values from `dsp_qdrift_flt_optimization` are passed to the flat-top fit of LegendSpecFits.jl, which uses them to correct for incomplete charge collection before the width is evaluated.

### Quality-Control Table

`dsp_qc_flt_optimization` and `dsp_qc_flt_optimization_compressed` run a short chain on the waveforms, typically a sample of events in a calibration line: baseline, PZ correction, `t50` from half of the PZ-corrected maximum, and the energy from the trapezoid with the default `rt` and `ft` of `config.default_flt_param.trap`. They return a table with the columns `energy`, `blmean`, `blslope`, `t50` and `qc_label`. A selection on these columns (for example on the baseline) is applied before the scans. The fourth argument is the classifier function of [Quality-Cut Label](@ref) or `missing`; with `missing` the label is `-1` and no classifier is evaluated. The two variants differ only in the preprocessing that the classifier receives (5 and 2 Haar passes).

### Q-drift

`dsp_qdrift_flt_optimization(wvfs, blmean, config, τ)` subtracts the given baseline means (`blmean`, one per waveform, for example the `blmean` column of the quality-control table), corrects the PZ, finds t0 and returns the Q-drift parameter, with the same definition as the column `qdrift` of [`dsp_icpc`](@ref). The values depend on `t0_threshold`, `qdrift_int_length` and `kwargs_pars`. In the compressed data layout the function is called with the windowed waveforms and the baseline means divided by `presum_rate`.

## Prepared Workflow

The scans above apply a filter at every grid point. When the same waveforms are evaluated many times, for example to search the (`rt`, `ft`) plane with an optimizer instead of a grid, the preparation can be done once:

1. [`dsp_flt_optimization_prepare`](@ref)`(wvfs, config, τ)` returns a [`FltOptPrep`](@ref) with the baseline-subtracted, PZ-corrected waveforms (`wvfs_pz`), `t50`, the baseline statistics (`bl_stats`), the signal estimator, and the sampling period.
2. [`dsp_energy_rtft`](@ref)`(prep, config, rt, ft, flt)` returns the energy of every waveform at one rise time and flat-top time. `flt` is `Val(:trap)`, `Val(:cusp)`, `Val(:zac)` or the corresponding `Symbol`. The energy is read at `t50 + rt + ft/2` for the trapezoid and at `t50 + flt_length/2` for CUSP and ZAC, as in the flat-top scans; `flt_length` is `config.flt_length_cusp` or `config.flt_length_zac`. `rt` and `ft` can be given in any time unit and number type. They are converted to floating point of the type of the sampling period (trapezoid) or of the filter length (CUSP, ZAC), so CUSP and ZAC accept filter times that are not multiples of the sampling period. The result is a `Vector{Float64}` for `Float64` waveforms.
3. The rise-time scans accept the prepared input as `dsp_{trap,cusp,zac}_rt_optimization(prep, config; ft)`. The flat-top scans are built on `dsp_energy_rtft` and prepare the waveforms once per call.

A `FltOptPrep` cannot be printed: displaying it raises a `StackOverflowError` in the display of its `SignalEstimator`. End assignments with a semicolon.

### Realized Flat-Top Time

The filters of RadiationDetectorDSP.jl work on whole samples, so the flat-top time that a filter applies can differ from the requested one. [`realized_flattop`](@ref)`(ft, flt_length, sampling_period, flt)` returns the applied value:

* Trapezoid: the flat top is `round(ft / sampling_period)` samples. `flt_length` is not used.
* CUSP and ZAC: the flat top is `round(ft / sampling_period)` samples, plus one sample if the filter length minus the flat top, both in samples, is odd.

The result is a multiple of the sampling period and differs from `ft` by less than two samples. Two requested flat-top times with the same realized value give the same filter and the same energies, so the realized value is the quantity to plot against when the energy resolution is evaluated as a function of the flat-top time.

### Cost and Threading

`dsp_energy_rtft` only reads `prep`, so it is safe to call from several threads on the same `prep`. Each call allocates about twice (trapezoid) or five times (CUSP, ZAC) the memory of `prep.wvfs_pz`. The number of simultaneous calls has to be bounded accordingly. The grid scans loop over the grid points one after the other.

## Current Filter Scans

`dsp_sg_optimization` and `dsp_mwa_optimization` compute A/E for every window length of the grid. They take the optimal trapezoid parameters from `pars_filter.trap.rt` and `pars_filter.trap.ft`; no default is used here, so both must be present. The energy is the trapezoid energy at `t50 + rt + ft/2`. The current is the maximum in `current_window` of

* the Savitzky-Golay derivative (degree `config.sg_flt_degree`, window `wl`) for `dsp_sg_optimization`, or
* the simple derivative smoothed by `config.a_mwa_passes` passes of a moving average of length `wl` for `dsp_mwa_optimization`.

The returned table has the columns

| Column | Content |
|:-------|:--------|
| `aoe` | For every waveform a vector with one A/E value per window length of the grid (a `VectorOfSimilarVectors`) |
| `energy` | Trapezoid energy |
| `blmean`, `blslope` | Baseline mean and slope |
| `t50` | Time of half the PZ-corrected maximum |
| `qc_label` | Classifier label, `-1` without a classifier |

The compressed variants `dsp_sg_optimization_compressed(wvfs_wdw, wvfs_pre, config, τ, pars_filter; presum_rate, f_evaluate_qc)` take the windowed waveforms first and the presummed waveforms second. The baseline, the energy and `t50` come from the presummed waveforms, and the current comes from the windowed waveforms. The keyword `presum_rate` is the presum rate of the data (the default is 8, so pass the actual value). The baseline mean of the presummed waveform divided by `presum_rate` is subtracted from the windowed waveform. The A/E values of the compressed variants therefore contain the energy in units of the presummed waveform.

The pulse-shape separation of events with a single energy deposit from events with multiple deposits is evaluated on A/E values of two calibration lines; see the [A/E Pulse-Shape Parameter](@ref) tutorial.

## A Synthetic Example

The example creates waveforms with a constant amplitude and a 1 µs linear rise, white noise, and an exponential tail, builds a [`DSPConfig`](@ref), and runs a rise-time scan and a flat-top scan of the trapezoidal filter. The grids are coarser than typical production grids to keep the example fast.

```@example filter_optimization
using LegendDSP
using PropDicts, Unitful, StableRNGs, Statistics, RadiationDetectorSignals, ArraysOfArrays, TypedTables

grid(start, stop, step) = PropDict(:start => start, :stop => stop, :step => step)
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
    :e_grid_trap => PropDict(:rt => grid(1.0u"µs", 13.0u"µs", 2.0u"µs"), :ft => grid(0.5u"µs", 4.5u"µs", 1.0u"µs")),
    :e_grid_zac => PropDict(:rt => grid(1.0u"µs", 13.0u"µs", 2.0u"µs"), :ft => grid(0.5u"µs", 4.5u"µs", 1.0u"µs")),
    :e_grid_cusp => PropDict(:rt => grid(1.0u"µs", 13.0u"µs", 2.0u"µs"), :ft => grid(0.5u"µs", 4.5u"µs", 1.0u"µs")),
    :a_grid_wl_sg => grid(30.0u"ns", 350.0u"ns", 64.0u"ns"),
    :a_grid_wl_mwa => grid(32.0u"ns", 320.0u"ns", 96.0u"ns"),
    :a_mwa_windows => [48.0u"ns"], :a_mwa_passes => 3, :a_sg_windows => [60.0u"ns"],
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

rng = StableRNG(3)
τ = 450.0u"µs"

function synthetic_waveform(amplitude; n = 8192, dt = 16.0u"ns", t_rise = 50.0u"µs", rise = 1.0u"µs", σ = 3.0)
    t = range(0.0u"ns", step = dt, length = n)
    x = ustrip.(u"µs", t) .- ustrip(u"µs", t_rise)
    shape = map(x) do xi
        xi < 0 ? 0.0 : xi < ustrip(u"µs", rise) ? xi / ustrip(u"µs", rise) :
            exp(-(xi - ustrip(u"µs", rise)) / ustrip(u"µs", τ))
    end
    RDWaveform(t, 1000 .+ amplitude .* shape .+ σ .* randn(rng, n))
end

wvfs = ArrayOfRDWaveforms([synthetic_waveform(5000.0) for _ in 1:60])
nothing # hide
```

The rise-time scan returns a matrix with one row per rise time of `config.e_grid_rt_trap`. The standard deviation of each row is the noise of the filter output at that rise time. For white noise it decreases with the rise time:

```@example filter_optimization
enc_grid = dsp_trap_rt_optimization(wvfs, config, τ)
size(enc_grid), collect(config.e_grid_rt_trap)
```

```@example filter_optimization
Table(rt = collect(config.e_grid_rt_trap), noise = vec(std(enc_grid; dims = 2)))
```

The flat-top scan at a chosen rise time returns the energy of every waveform for every flat-top time. All waveforms have the same amplitude here, so the standard deviation of a row is the resolution contribution of the noise at that flat-top time:

```@example filter_optimization
e_grid = dsp_trap_ft_optimization(wvfs, config, τ, 9.0u"µs")
Table(ft = collect(config.e_grid_ft_trap), mean_energy = vec(mean(e_grid; dims = 2)), width = vec(std(e_grid; dims = 2)))
```

The prepared workflow does the decoding, baseline subtraction and PZ correction once. [`dsp_energy_rtft`](@ref) then evaluates single points, and the rise-time scan runs on the prepared input:

```@example filter_optimization
prep = dsp_flt_optimization_prepare(wvfs, config, τ);
enc_grid_prep = dsp_trap_rt_optimization(prep, config)
enc_grid_prep == enc_grid
```

The energy at a single (`rt`, `ft`) point equals the corresponding row of the flat-top scan, apart from the conversion of the scan result to `Float32`:

```@example filter_optimization
ft_grid = collect(config.e_grid_ft_trap)
e_point = dsp_energy_rtft(prep, config, 9.0u"µs", ft_grid[2], :trap)
(eltype(e_point), maximum(abs.(e_point .- e_grid[2, :])) < 1e-3 * mean(e_point))
```

The filter times are not limited to the grid. Several points can be evaluated from different threads on the same `prep`:

```@example filter_optimization
points = [(rt, ft) for rt in (5.0u"µs", 9.0u"µs"), ft in (1.0u"µs", 2.0u"µs", 3.0u"µs")]
widths = zeros(size(points))
Threads.@threads for i in eachindex(points)
    widths[i] = std(dsp_energy_rtft(prep, config, points[i]..., :cusp))
end
Table(rt = first.(vec(points)), ft = last.(vec(points)), width = vec(widths))
```

The flat-top time that a filter applies is a multiple of the sampling period of 16 ns. For the trapezoid it is the rounded number of samples; for CUSP and ZAC it also depends on the filter length, because the filter length minus the flat top has to be an even number of samples:

```@example filter_optimization
dt = prep.sampling_period
Table(requested = ft_grid,
      trap = realized_flattop.(ft_grid, config.flt_length_cusp, dt, :trap),
      cusp = realized_flattop.(ft_grid, config.flt_length_cusp, dt, :cusp))
```

The quality-control table of the same waveforms, without a classifier:

```@example filter_optimization
qc = dsp_qc_flt_optimization(wvfs, config, τ, missing)
qc[1:3]
```

The A/E scan needs the trapezoid parameters in `pars_filter`. To get a spread of A/E, the waveforms below have different amplitudes. A/E for every window length is stored in the column `aoe`:

```@example filter_optimization
pars_filter = PropDict(:trap => PropDict(:rt => 9.0u"µs", :ft => 2.5u"µs"))
wvfs_amp = ArrayOfRDWaveforms([synthetic_waveform(a) for a in range(3000.0, 8000.0, length = 40)])

sg = dsp_sg_optimization(wvfs_amp, config, τ, pars_filter)
mwa = dsp_mwa_optimization(wvfs_amp, config, τ, pars_filter)

aoe_sg = reduce(hcat, sg.aoe)    # window lengths in rows, waveforms in columns
aoe_mwa = reduce(hcat, mwa.aoe)
relative_spread(aoe) = vec(std(aoe; dims = 2) ./ mean(aoe; dims = 2))

(sg = Table(wl = collect(config.a_grid_wl_sg), relative_spread = relative_spread(aoe_sg)),
 mwa = Table(wl = collect(config.a_grid_wl_mwa), relative_spread = relative_spread(aoe_mwa)))
```

All waveforms have the same rise time, so the relative spread of A/E is the noise contribution of the current filter, which decreases with the window length until the filter smooths the rise itself. For real data the scan compares the A/E distributions of two calibration lines instead.

## DPLMS Filter

### Method

A DPLMS filter is a finite-impulse-response energy filter whose weights are the solution of a constrained least-squares problem instead of a fixed shape such as a trapezoid. The weights are derived directly from measured signal and noise waveforms. The method is described by D'Andrea et al., [arXiv:2209.09863](https://arxiv.org/abs/2209.09863). The weights `x` minimize the filter output noise `xᵀ N x`, with the noise matrix `N` measured on baseline samples, under the constraints

* the response to the reference pulse `r` is one, `x ⋅ r = 1`,
* the output is flat around the read-out time (flat top), so that the energy does not depend on small shifts of the pulse, and
* the coefficients sum to zero (zero area), so that a constant baseline offset gives no output.

The constraints enter as quadratic penalties with weights, which turns the problem into one linear system. A further weight `zo` can make the area of the output for a step input zero, as for a ZAC filter. [`dplms_synthesis`](@ref) documents the system and the weights.

### Input

The functions work on presummed waveforms that are baseline subtracted and pole-zero corrected (the waveforms `wvfs_pre` inside [`dsp_icpc_compressed`](@ref)). The waveforms are aligned on `t50`, the time at which they reach 50 % of the maximum of a `TrapezoidalChargeFilter(10u"µs", 10u"µs")`; this is the column `t50_pre` of the compressed chain. The reference pulse is the mean of the aligned waveforms (superpulse), and its amplitude after the rise is the unit of all matrices. The noise matrix needs at least as many baseline samples inside `config.bl_window` as the filter has coefficients.

### Building Blocks

The functions of `dplms.jl` are the steps of the synthesis:

| Function | Result |
|:---------|:-------|
| [`align_waveform`](@ref), [`align_waveforms`](@ref) | Waveforms shifted in time by a fractional number of samples with a cubic interpolation, so that the feature at `t_align` (here `t50`) lands at a common time |
| [`aligned_windows`](@ref) | Matrix of the samples from `n_pre` before to `n_post` after the common time, one column per waveform |
| [`dplms_noise_matrix`](@ref) | Noise matrix: the mean outer product of all windows of the filter length inside the baseline window |
| [`dplms_reference`](@ref) | Superpulse of the windows, normalized to unit amplitude on a plateau after the rise |
| [`dplms_flattop_matrix`](@ref) | Matrix whose quadratic form is the mean squared slope of the filter output, from reference pulses widened by `ff` samples on each side |
| [`dplms_binned_reference`](@ref) | References per drift-time bin, as the matrix and vector of an additional constraint that the filter has unit response for every bin |
| [`dplms_synthesis`](@ref) | Solution of the linear system for the weights `nm`, `za`, `ft`, `ct`, `zo`: the filter, its predicted noise, and its area |

The filter returned by `dplms_synthesis` is the impulse response `coeffs`. The filter output at the last sample of a window of `n_pre + n_post + 1` samples is the energy, so the read-out time after `t50` is `n_post` samples.

### Optimization over a Grid

[`dsp_dplms_optimization`](@ref)`(wvfs, config, τ, drift; kwargs...)` runs all steps on a set of presummed waveforms, for example events of a calibration line, and synthesizes one filter for every point of a grid. `drift` is a drift-time proxy of every waveform; its quantiles define the `n_bins` drift-time bins of the reference constraint. The function subtracts the baseline in `config.bl_window`, corrects the pole-zero with `InvCRFilter(τ)`, determines `t50`, aligns the waveforms to the sample nearest the median `t50`, and builds the matrices. Unless given, the window has the longest symmetric extent `n_pre = n_post` that leaves `n_noise_windows` windows per event inside the baseline window for the noise matrix.

A grid point is a `NamedTuple` with the weights `nm` (noise), `za` (zero area), `ft` (flat top), `ct` (binned references), `zo` (output area) and the flat-top widening `ff` in samples. [`dplms_default_grid`](@ref)`(n_rise; kwargs...)` builds the grid with the widenings `ff = k·n_rise` for `k` in 1:4 and the weights given as keywords. `n_rise` is the number of samples from the alignment point until the slowest drift-time reference reaches 99.9 % of its plateau. The widening sets how far the pulse can move while the output stays flat; larger values give a longer flat top and more noise. Only the ratios of the weights matter. The keyword `grid` takes a vector of grid points or a function of `n_rise` that returns one, for example `grid = n -> dplms_default_grid(n; za = 0)` for filters without the zero-area condition.

The function returns a `NamedTuple`. The energies `e_grid` and filter-output maxima `emax_grid` have one row per grid point and one column per waveform. The filters are in `coeffs` with the common read-out time `t_pickoff`. For every grid point `noise_out` is the predicted output noise in units of the reference amplitude, `zero_area` the coefficient sum, `output_area` the area of the output for a step input, and `flat_top` and `shaping_time` the length of the flat top and of the rise of the output for the drift-time references. The matrices, references, superpulses and the response of every filter to the superpulse are returned as well; see the docstring for the list. As for the other scans, choosing the best grid point from the widths of the energy distributions is a separate step.

### Energy in the Compressed Chain

[`dsp_dplms_energy`](@ref) applies a filter to baseline-subtracted, pole-zero-corrected waveforms and reads the value at `t50 + t_pickoff` with a signal estimator. [`dsp_icpc_compressed`](@ref) calls it on the presummed waveforms with the filter from [`get_dplms_pars`](@ref)`(pars_filter, config, Δt)` and the column `t50_pre`. It stores the energy as `e_dplms` and the maximum of the output in `config.current_window` shifted by `t_pickoff` as `e_dplms_max`.

`get_dplms_pars` returns `coefficients` and `t_pickoff` of `pars_filter.dplms` if the entry exists. Otherwise it builds the trapezoid of `config.default_flt_param.dplms` (`rt` and `ft`) as a filter at the sampling interval `Δt`, read out at `rt + ft/2`, so that `e_dplms` is the energy of `TrapezoidalChargeFilter(rt, ft)`. Without `pars_filter.dplms` and without `default_flt_param.dplms` the function throws an `ArgumentError`. The production dataflow does not call the DPLMS synthesis, see [Use in Juleana.jl](@ref).

### A Synthetic DPLMS Example

The example synthesizes one DPLMS filter from synthetic presummed waveforms. The pulses are steps with rise times between 0.5 and 2.5 µs at a jittering onset, with white noise and a random walk as baseline noise. There is no decay tail, so no pole-zero correction is needed. The rise time stands in for the drift-time proxy. The filter is compared with a trapezoid of the same length on the same events. The example only illustrates the functions; the noise model is not that of a detector.

```@example dplms
using LegendDSP
using RadiationDetectorDSP, RadiationDetectorSignals, ArraysOfArrays
using Unitful, StableRNGs, Statistics, IntervalSets, TypedTables

rng = StableRNG(7)
dt = 0.5u"µs"
t = range(0.0u"µs", step = dt, length = 600)

function synthetic_pulse(rng, t; amplitude = 1000.0, σ_white = 3.0, σ_walk = 0.6)
    t0 = 100.0 + 3 * randn(rng)
    rise = 0.5 + 2 * rand(rng)
    shape = clamp.((ustrip.(u"µs", t) .- t0) ./ rise, 0, 1)
    noise = σ_white .* randn(rng, length(t)) .+ σ_walk .* cumsum(randn(rng, length(t)))
    (wvf = RDWaveform(t, 20.0 .+ amplitude .* shape .+ noise), rise = rise)
end

pulses = [synthetic_pulse(rng, t) for _ in 1:400]
drift = [p.rise for p in pulses]
bl_window = 0.0u"µs" .. 80.0u"µs"
wvfs = ArrayOfRDWaveforms([p.wvf for p in pulses])
wvfs = shift_waveform.(wvfs, -[mean(w.signal[1:161]) for w in wvfs])   # baseline subtraction
nothing # hide
```

`t50` is the time of half the maximum of the 10/10 µs trapezoid. The waveforms are aligned to the sample nearest the median `t50`, and a window of 101 samples (50 µs) around it gives the reference and the noise matrix. The flat-top matrix uses a widening of 6 samples (3 µs), and the drift-time bins are terciles of the rise time:

```@example dplms
e_101010 = extremestats.(TrapezoidalChargeFilter(10u"µs", 10u"µs").(wvfs)).max
t50 = get_threshold(wvfs, 0.5 .* e_101010; mintot = 2.0u"µs")

n_pre, n_post, ff = 40, 60, 6
n = n_pre + n_post + 1
t_target = t[findmin(abs.(t .- median(t50)))[2]]
wvfs_aligned = align_waveforms(wvfs, t50, t_target)

windows = aligned_windows(wvfs_aligned, t_target, n_pre, n_post)
ref, plateau = dplms_reference(windows, (n - 19):n)
noise = dplms_noise_matrix(wvfs, bl_window, n) ./ plateau^2
bins = cld.(invperm(sortperm(drift)) .* 3, length(drift))
binned = dplms_binned_reference(windows, bins, plateau)
wide = aligned_windows(wvfs_aligned, t_target, n_pre + ff, n_post + ff)
flattop = dplms_flattop_matrix(dplms_binned_reference(wide, bins, plateau).refs, ff)

flt = dplms_synthesis(noise, ref, flattop, binned; nm = 1.0, za = 1e3, ft = 1.0, zo = 1e3)
(noise_out = flt.noise_out, zero_area = flt.zero_area, output_area = flt.output_area)
```

The predicted output noise is in units of the plateau amplitude. The coefficient sum and the area of the output for a step are zero up to rounding. The energies of the DPLMS filter and of a trapezoid with 20 µs rise time and 10 µs flat top are read out with [`dsp_dplms_energy`](@ref). The trapezoid is given by the coefficients that [`get_dplms_pars`](@ref) builds for the default filter:

```@example dplms
estimator = SignalEstimator(PolynomialDNI(3, 2.0u"µs"))
current_window = 80.0u"µs" .. 120.0u"µs"
n_rt, n_ft = 40, 20
trap = [fill(1 / n_rt, n_rt); zeros(n_ft); fill(-1 / n_rt, n_rt)]

e_dplms = dsp_dplms_energy(wvfs, flt.coeffs, t50, n_post * dt, current_window, estimator).e
e_trap = dsp_dplms_energy(wvfs, trap, t50, (n_rt + n_ft / 2) * dt, current_window, estimator).e

Table(filter = ["DPLMS", "trapezoid"], mean = [mean(e_dplms), mean(e_trap)], width = [std(e_dplms), std(e_trap)])
```

All pulses have the same amplitude, so the width is the noise of the energy estimate. The width of the DPLMS energies agrees with `noise_out` times the plateau amplitude. For this noise, which has a random-walk component, the DPLMS filter gives the smaller width. The weights show the difference: the trapezoid weights the samples of its two averaging windows equally, whereas the DPLMS weights follow the noise and the pulse shape.

```@example dplms
using CairoMakie
CairoMakie.activate!(type = "svg")

# time of each weight relative to t50 when the filter is read out at t50 + t_pickoff
t_weights(c, t_pickoff) = ustrip.(u"µs", t_pickoff .- (0:length(c)-1) .* dt)

fig = Figure(size = (700, 320))
ax = Axis(fig[1, 1], xlabel = "time of the weighted sample relative to t50 (µs)", ylabel = "weight")
lines!(ax, t_weights(flt.coeffs, n_post * dt), flt.coeffs, label = "DPLMS")
lines!(ax, t_weights(trap, (n_rt + n_ft / 2) * dt), trap, label = "trapezoid")
axislegend(ax, position = :rt)
fig
```

## Role in the Production Dataflow

In the production dataflow the scans form two steps that run per detector on waveforms of a calibration gamma line, after a selection with the quality-control table:

1. **Energy filter.** For each filter in the configuration (`trap` and `cusp` in the default configurations), the rise-time scan runs with a fixed flat-top time from the configuration. The width of the filter-output noise as a function of the rise time gives the optimal rise time. The flat-top scan at this rise time, with the Q-drift values from `dsp_qdrift_flt_optimization`, gives the energy resolution as a function of the flat-top time and so the optimal flat-top time.
2. **Current filter.** The scans `dsp_sg_optimization_compressed` and `dsp_mwa_optimization_compressed` run on waveforms of a double-escape and a single-escape peak. The window length is chosen from the A/E distributions of the two peaks by `fit_sf_wl` of LegendSpecFits.jl and stored.

The dataflow calls the rise-time and flat-top scans with the waveforms and `τ`; it does not use the prepared workflow. The results are written to the `fltopt` and `aoeopt` parameter sets and merged into the `pars_filter` argument of `dsp_icpc_compressed`. The function names are built from the configuration keys; this contract is described in [Use in Juleana.jl](@ref).
