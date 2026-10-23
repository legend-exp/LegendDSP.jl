# Filter Optimization

The energy filters (trapezoidal, CUSP, ZAC) and the current filters (Savitzky-Golay, moving window average) of the HPGe chains have free parameters: the rise time `rt` and flat-top time `ft` of an energy filter and the window length `wl` of a current filter. The optimal values depend on the detector and its electronics. LegendDSP.jl provides functions that evaluate a filter on a grid of parameter values for a set of waveforms of one detector. They return the filter output for every grid point and every waveform. The statistical evaluation of these grids (noise width, energy resolution, A/E separation) is done by the fit functions of LegendSpecFits.jl; see the [LegendSpecFits.jl documentation](https://legend-exp.github.io/LegendSpecFits.jl/stable/). The optimal values enter the DSP chains as the argument `pars_filter` ([HPGe DSP Chains](@ref)).

All scans start with the same preparation: encoded waveforms are decoded, the baseline mean in `config.bl_window` is subtracted, and the waveforms are corrected for the pole-zero with `InvCRFilter(τ)` (see [Decay Time](@ref)). The scans therefore need the decay time `τ` of the detector. The grids and the other settings come from the [`DSPConfig`](@ref).

## Overview

| Function | Scanned parameter | Grid in `config` | Result |
|:---------|:------------------|:-----------------|:-------|
| [`dsp_trap_rt_optimization`](@ref), [`dsp_cusp_rt_optimization`](@ref), [`dsp_zac_rt_optimization`](@ref) | rise time `rt` at fixed flat-top time (keyword `ft`, default 2 µs) | `e_grid_rt_trap`, `e_grid_rt_cusp`, `e_grid_rt_zac` | `Matrix{Float64}`, rows: grid points, columns: waveforms |
| [`dsp_trap_ft_optimization`](@ref), [`dsp_cusp_ft_optimization`](@ref), [`dsp_zac_ft_optimization`](@ref) | flat-top time `ft` at fixed rise time `rt` (argument) | `e_grid_ft_trap`, `e_grid_ft_cusp`, `e_grid_ft_zac` | `Matrix{Union{Missing,Float32}}`, rows: grid points, columns: waveforms |
| [`dsp_qc_flt_optimization`](@ref), [`dsp_qc_flt_optimization_compressed`](@ref) | none | none | `Table` with the energy at the default trapezoid, baseline mean and slope, `t50`, and `qc_label` |
| [`dsp_qdrift_flt_optimization`](@ref) | none | `qdrift_int_length` | Vector of Q-drift values |
| [`dsp_sg_optimization`](@ref), [`dsp_sg_optimization_compressed`](@ref) | window length of the Savitzky-Golay current filter | `a_grid_wl_sg` | `Table` with A/E for every window length |
| [`dsp_mwa_optimization`](@ref), [`dsp_mwa_optimization_compressed`](@ref) | window length of the moving-window-average current filter | `a_grid_wl_mwa` | `Table` with A/E for every window length |

There are no `_compressed` variants of the rise-time and flat-top scans. They work on a single waveform per event, and the dataflow passes them the presummed waveforms.

## Energy Filter Scans

### Rise Time

`dsp_{trap,cusp,zac}_rt_optimization(wvfs, config, τ; ft)` applies the filter with every rise time in the grid to all waveforms and reads the filter output at the fixed time `config.enc_pickoff_{trap,cusp,zac}`. This pick-off time lies in the part of the trace before the signal rise, so the values are samples of the filter-output noise and no signal. Their spread over the waveforms, as a function of the rise time, is the electronic noise (ENC) of the filter. The function returns the values; the width is determined by the fit. For the CUSP and ZAC filters the grid value is the shaping time σ, and `ft` is the length of the flat top; the total filter length is `config.flt_length_cusp` or `config.flt_length_zac`.

### Flat-Top Time

`dsp_{trap,cusp,zac}_ft_optimization(wvfs, config, τ, rt)` applies the filter with the given rise time `rt` and every flat-top time in the grid. The energy is read at a pick-off time that depends on the filter: `t50 + rt + ft/2` for the trapezoid and `t50 + flt_length/2` for CUSP and ZAC, where `t50` is the time at which the PZ-corrected waveform reaches half of its maximum. The result is the energy of every waveform for every flat-top time. A flat-top time shorter than the rise time of the signal gives a lower energy (ballistic deficit). For a calibration line with many events, the width of the line as a function of the flat-top time shows the flat-top time that gives the best resolution. The Q-drift values from `dsp_qdrift_flt_optimization` are passed to the flat-top fit of LegendSpecFits.jl, which uses them to correct for incomplete charge collection before the width is evaluated.

### Quality-Control Table

`dsp_qc_flt_optimization` and `dsp_qc_flt_optimization_compressed` run a short chain on the waveforms, typically a sample of events in a calibration line: baseline, PZ correction, `t50` from half of the PZ-corrected maximum, and the energy from the trapezoid with the default `rt` and `ft` of `config.default_flt_param.trap`. They return a table with the columns `energy`, `blmean`, `blslope`, `t50` and `qc_label`. A selection on these columns (for example on the baseline) is applied before the scans. The fourth argument is the classifier function of [Quality-Cut Label](@ref) or `missing`; with `missing` the label is `-1` and no classifier is evaluated. The two variants differ only in the preprocessing that the classifier receives (5 and 2 Haar passes).

### Q-drift

`dsp_qdrift_flt_optimization(wvfs, blmean, config, τ)` subtracts the given baseline means (`blmean`, one per waveform, for example the `blmean` column of the quality-control table), corrects the PZ, finds t0 and returns the Q-drift parameter, with the same definition as the column `qdrift` of [`dsp_icpc`](@ref). The values depend on `t0_threshold`, `qdrift_int_length` and `kwargs_pars`. In the compressed data layout the function is called with the windowed waveforms and the baseline means divided by `presum_rate`.

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

## Role in the Production Dataflow

In the production dataflow the scans form two steps that run per detector on waveforms of a calibration gamma line, after a selection with the quality-control table:

1. **Energy filter.** For each filter in the configuration (`trap` and `cusp` in the default configurations), the rise-time scan runs with a fixed flat-top time from the configuration. The width of the filter-output noise as a function of the rise time gives the optimal rise time. The flat-top scan at this rise time, with the Q-drift values from `dsp_qdrift_flt_optimization`, gives the energy resolution as a function of the flat-top time and so the optimal flat-top time.
2. **Current filter.** The scans `dsp_sg_optimization_compressed` and `dsp_mwa_optimization_compressed` run on waveforms of a double-escape and a single-escape peak. The window length is chosen from the A/E distributions of the two peaks by `fit_sf_wl` of LegendSpecFits.jl and stored.

The results are written to the `fltopt` and `aoeopt` parameter sets and merged into the `pars_filter` argument of `dsp_icpc_compressed`. The function names are built from the configuration keys; this contract is described in [Use in Juleana.jl](@ref).
