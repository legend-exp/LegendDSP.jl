# Waveforms and Configuration

All processing functions of LegendDSP.jl take digitized detector waveforms and a set of configuration parameters. This page describes the waveform types, the two layouts in which waveforms arrive from the raw data tier (uncompressed and compressed), and the [`DSPConfig`](@ref) type that carries the parameters of the HPGe processing.

## Waveform Types

A single waveform is a `RDWaveform` from [RadiationDetectorSignals.jl](https://github.com/JuliaPhysics/RadiationDetectorSignals.jl). It pairs a time axis with a vector of samples. The time axis is a range of Unitful quantities, so that all window lengths, filter lengths and time points in LegendDSP.jl are quantities as well (`µs`, `ns`). The samples are plain numbers in ADC units.

```@example waveforms_config
using LegendDSP
using RadiationDetectorSignals
using Unitful, StableRNGs

rng = StableRNG(42)
t = range(0.0u"ns", step = 16.0u"ns", length = 4096)
signal = 1000 .+ 3 .* randn(rng, length(t))
wf = RDWaveform(t, signal)
(length(wf.signal), step(wf.time))
```

A set of waveforms is stored as an `ArrayOfRDWaveforms`. It keeps all samples in one matrix, all time axes in one vector, and behaves like a vector of `RDWaveform`. This is the type of the `waveform` columns that the DSP chains read:

```@example waveforms_config
using ArraysOfArrays

wvfs = ArrayOfRDWaveforms([RDWaveform(t, 1000 .+ 3 .* randn(rng, length(t))) for _ in 1:4])
(length(wvfs), eltype(wvfs) <: RDWaveform, step(wvfs[1].time))
```

The DSP functions broadcast over waveform sets. Almost all of them treat each waveform independently and return one value (or one small vector) per waveform. Time points returned by LegendDSP.jl are `Unitful.Time` quantities in the unit of the time axis or in µs, as listed for each function in [HPGe DSP Chains](@ref) and [Waveform Parameters](@ref).

Functions that use a window on the time axis (`bl_window`, `tail_window`, `current_window`) convert the time interval into sample indices with the step of the time axis of the waveform they are applied to. Some functions (`get_triggers`, the filters scaled by `flt_length_cusp`, `flt_length_zac`) read the step and the first time of the first waveform of a set and apply it to all waveforms. All waveforms of one call must therefore share the same time axis.

## Uncompressed and Compressed Input

The HPGe chains exist in two variants that differ in the layout of the raw waveform columns they read.

| Variant | Waveform columns | Further columns |
|:--------|:-----------------|:----------------|
| [`dsp_icpc`](@ref), [`dsp_puls`](@ref), [`dsp_aux`](@ref) | `waveform` | `baseline`, `timestamp`, `eventnumber`, `daqenergy` |
| [`dsp_icpc_compressed`](@ref), [`dsp_puls_compressed`](@ref), [`dsp_aux_compressed`](@ref) | `waveform_presummed`, `waveform_windowed` (the pulser and auxiliary chains read only `waveform_presummed`) | `presum_rate`, `baseline`, `timestamp`, `eventnumber`, `daqenergy`, `t_sat_lo`, `t_sat_hi`, `deadtime` (the last three only for `dsp_icpc_compressed`) |

The compressed layout stores each event as two waveforms:

* `waveform_presummed` covers the whole trace. Each sample is the sum of `presum_rate` consecutive ADC samples, so the time step is `presum_rate` times the sampling period and the amplitude is `presum_rate` times the amplitude of the original samples. It is used for baseline, tail, energy filters and the saturation counts.
* `waveform_windowed` covers a short interval around the signal rise at the full sampling rate and the original amplitude scale. It is used for timing points (t0, t10 to t99), the current (A), the Q-drift and the LQ parameters. Its time axis must use the same origin as the time axis of `waveform_presummed`: time points found on one waveform are applied to the other, and the windows of the configuration are interpreted on both.

All `_compressed` chains subtract the baseline mean measured on the presummed waveform from both waveforms; the windowed waveform is shifted by the mean divided by `presum_rate`. Energies derived from the presummed waveform are therefore in units of `presum_rate` times ADC counts, and energies derived from the windowed waveform are in ADC counts. The chains read `presum_rate` as a single value per call: all events of one table must have the same presum rate.

Waveform columns may be stored encoded (`EncodedArrays.jl` with a lossless codec). Every DSP function in LegendDSP.jl calls `decode_data` from LegendDataTypes.jl on its input first. `decode_data` decodes encoded arrays inside tables, named tuples, `RDWaveform` and `ArrayOfRDWaveforms`, and returns any other input unchanged. Both the plain and the `_compressed` functions therefore accept encoded and already decoded columns. The `_compressed` suffix selects the column layout, not the decoding.

The SiPM variant `dsp_sipm_compressed` reads `waveform_bit_drop` instead; see [SiPM and PMT Processing](@ref).

## The Configuration Type DSPConfig

The HPGe, pulser and auxiliary chains, the decay-time extraction and the filter optimization scans read their settings from a [`DSPConfig`](@ref). It is an immutable struct with 28 fields. All lengths are Unitful quantities, all grids are ranges of quantities.

| Field | Content | Read by |
|:------|:--------|:--------|
| `enc_pickoff_trap`, `enc_pickoff_zac`, `enc_pickoff_cusp` | Time at which the filtered waveform is sampled for the noise (ENC) scan of the rise time | `dsp_{trap,zac,cusp}_rt_optimization` |
| `flt_length_cusp`, `flt_length_zac` | Total length of the CUSP and ZAC filter kernels | `dsp_icpc(_compressed)`, `dsp_{cusp,zac}_{rt,ft}_optimization` |
| `t0_threshold` | Threshold (ADC counts) on the fast asymmetric trapezoid for the t0 determination | `dsp_icpc(_compressed)`, `dsp_qdrift_flt_optimization` |
| `inTraceCut_std_threshold` | In-trace pile-up threshold in multiples of the baseline standard deviation | `dsp_icpc(_compressed)` |
| `sg_flt_degree` | Polynomial degree of the Savitzky-Golay current filter | `dsp_icpc(_compressed)`, `dsp_{sg,mwa}_optimization(_compressed)` |
| `bl_window` | Interval for the baseline fit | all chains and scans |
| `tail_window` | Interval for the tail fit (decay time) and the tail statistics after the pole-zero correction | `dsp_icpc(_compressed)`, `dsp_decay_times` |
| `current_window` | Interval in which the maximum current (A) is searched | `dsp_icpc(_compressed)`, `dsp_{sg,mwa}_optimization(_compressed)` |
| `qdrift_int_length` | Range whose first and last element delimit the Q-drift integration window (relative to t0) | `dsp_icpc(_compressed)`, `dsp_qdrift_flt_optimization` |
| `lq_windows` | List of `[early, gap, late]` triples, one set of LQ columns per triple | `dsp_icpc(_compressed)` |
| `e_grid_rt_trap`, `e_grid_ft_trap`, `e_grid_rt_zac`, `e_grid_ft_zac`, `e_grid_rt_cusp`, `e_grid_ft_cusp` | Rise-time and flat-top-time grids of the energy filter scans | `dsp_{trap,zac,cusp}_{rt,ft}_optimization` |
| `a_grid_wl_sg`, `a_grid_wl_mwa` | Window-length grids of the SG and MWA current filters | `dsp_sg_optimization(_compressed)`, `dsp_mwa_optimization(_compressed)` |
| `a_mwa_windows`, `a_sg_windows` | Fixed window lengths for which extra current columns `a_mwa_<ns>` and `a_sg_<ns>` are written | `dsp_icpc(_compressed)` |
| `a_mwa_passes` | Number of alternating moving-average passes of the MWA current filter | `dsp_icpc(_compressed)`, `dsp_mwa_optimization(_compressed)` |
| `default_flt_param` | Filter parameters used when a detector has no optimized value (see below) | `dsp_icpc(_compressed)`, `dsp_qc_flt_optimization(_compressed)` |
| `kwargs_pars` | Further settings, see below | most functions |
| `auxpz1_window`, `auxpz2_window` | Two additional intervals on the pole-zero corrected presummed waveform | `dsp_icpc_compressed` |

`dsp_puls`, `dsp_puls_compressed` use only `bl_window`. `dsp_aux` and `dsp_aux_compressed` accept a `DSPConfig` for a uniform call signature but do not read it.

`default_flt_param` is a `PropDict` with the entries `sg` and `mwa` (window lengths) and `trap`, `zac`, `cusp` (each with `rt` and `ft`). The chains look up optimized values in the `pars_filter` argument first and use `default_flt_param` for every entry that `pars_filter` lacks.

`kwargs_pars` is a `PropDict` with the following entries:

| Key | Content |
|:----|:--------|
| `fc_bit_depth` | Bit depth of the FADC; sets the upper saturation level |
| `t0_flt_pars` | Pre-rise averaging time, gap time and post-rise averaging time of the fast asymmetric trapezoid used for t0 (a vector of three times) |
| `t0_mintot` | Minimum time over threshold for the t0 pick-off |
| `tx_mintot` | Minimum time over threshold for the pick-offs t10 to t99 and `t50_current` |
| `intrace_mintot` | Minimum time over threshold for the in-trace pile-up search |
| `sig_interpolation_order`, `sig_interpolation_length` | Polynomial order and length of the signal estimator that reads the filter output at the pick-off time |
| `int_interpolation_order`, `int_interpolation_length` | Same, for the integrated waveform (Q-drift, LQ) |

### Building a DSPConfig from a PropDict

`DSPConfig(pd)` takes a `PropDict` with one key per setting. Quantities must already carry their units. Windows are sub-dictionaries with `min` and `max`; each grid is a sub-dictionary with `start`, `stop` and `step` (for `e_grid_trap`, `e_grid_zac`, `e_grid_cusp` one such triple for `rt` and one for `ft`). The constructor turns the triples into ranges and the windows into closed intervals. The following keys of `pd` differ from the field names:

* `flt_defaults` becomes the field `default_flt_param`.
* `qdrift_int_length` is a range or a two-element vector; the field becomes a range from the first to the last element in steps of one tenth of the unit.
* `a_grid_wl_mwa` is optional and falls back to `a_grid_wl_sg`.
* Each entry of `lq_windows` must have exactly three elements, otherwise an `ArgumentError` is thrown.

A synthetic example with generic values:

```@example waveforms_config
using PropDicts, IntervalSets

grid(start, stop, step) = PropDict(:start => start, :stop => stop, :step => step)
rt_ft_grid = PropDict(:rt => grid(1.0u"µs", 16.0u"µs", 0.5u"µs"),
                      :ft => grid(1.0u"µs", 4.0u"µs", 0.2u"µs"))

pd = PropDict(
    :enc_pickoff_trap => 40.0u"µs", :enc_pickoff_zac => 41.0u"µs", :enc_pickoff_cusp => 41.0u"µs",
    :flt_length_cusp => 38.0u"µs", :flt_length_zac => 38.0u"µs",
    :bl_window => PropDict(:min => 0.0u"µs", :max => 39.0u"µs"),
    :tail_window => PropDict(:min => 70.0u"µs", :max => 110.0u"µs"),
    :current_window => PropDict(:min => 43.0u"µs", :max => 62.0u"µs"),
    :auxpz1_window => PropDict(:min => 70.0u"µs", :max => 90.0u"µs"),
    :auxpz2_window => PropDict(:min => 90.0u"µs", :max => 110.0u"µs"),
    :t0_threshold => 4.0,
    :inTraceCut_std_threshold => 5,
    :sg_flt_degree => 3,
    :qdrift_int_length => [2.5u"µs", 5.0u"µs"],
    :lq_windows => [[2.5, 0.0, 2.5]u"µs"],
    :e_grid_trap => rt_ft_grid, :e_grid_zac => rt_ft_grid, :e_grid_cusp => rt_ft_grid,
    :a_grid_wl_sg => grid(30.0u"ns", 350.0u"ns", 32.0u"ns"),
    :a_mwa_windows => [48.0u"ns", 288.0u"ns"],
    :a_mwa_passes => 3,
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
)

config = DSPConfig(pd)
```

The fields hold typed values that the algorithms use directly:

```@example waveforms_config
config.bl_window, config.e_grid_rt_trap, config.a_grid_wl_mwa, config.qdrift_int_length
```

The same config object is used for every channel of one call. In the production dataflow, a per-detector config is the merge of a `default` section with a detector-specific section of the same layout; see [Use in Juleana.jl](@ref).

SiPM and PMT processing does not use `DSPConfig`; [`dsp_sipm`](@ref) and `dsp_pmts` take a plain `PropDict` whose keys are listed in [SiPM and PMT Processing](@ref).

### Optimized Filter Parameters

The functions [`dsp_icpc`](@ref) and [`dsp_icpc_compressed`](@ref) take a second `PropDict`, `pars_filter`, with the per-detector optimal filter parameters found by the scans described in [Filter Optimization](@ref). Its layout is

```
trap: {rt, ft}     cusp: {rt, ft}     zac: {rt, ft}     sg: {wl}     mwa: {wl}
```

with Unitful quantities as leaves. Entries may be missing; the corresponding value of `config.default_flt_param` is used instead. An empty `PropDict()` therefore runs the chain entirely with the defaults of the configuration.

```@example waveforms_config
pars_filter = PropDict(:trap => PropDict(:rt => 8.0u"µs", :ft => 2.0u"µs"))
config.default_flt_param.cusp.rt, haskey(pars_filter, :cusp)
```
