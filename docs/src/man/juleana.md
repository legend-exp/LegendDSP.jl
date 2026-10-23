# Use in Juleana.jl

[Juleana.jl](https://github.com/legend-exp/Juleana.jl) is the production dataflow of the LEGEND Julia analysis. Its DSP processors read raw waveforms of one channel at a time, call the functions of LegendDSP.jl, and store the returned tables in the DSP tier. Other processors determine the filter parameters and the decay time that the DSP processors need. This page describes how the dataflow builds the configuration, which processor calls which function, and how it selects LegendDSP.jl functions by name. The code examples cannot run without access to LEGEND data and metadata and are not executed.

## Building the Configuration

The DSP settings are part of the LEGEND production-data metadata. They are read through LegendDataManagement.jl for the validity period of a file key:

```julia
using LegendDataManagement, LegendDSP, PropDicts

l200 = LegendData(:l200)
filekey = start_filekey(l200, (:p03, :r000, :cal))

dsp_config_pd = dataprod_config(l200).dsp(filekey)   # PropDict
```

The metadata directory selects and merges the configuration files that are valid at the time of the file key. The result has a `default` section and optionally one section per detector, all with the same layout. Values with units are stored as `{val, unit}` pairs and are read as Unitful quantities. The DSP processors build the configuration of one detector by merging the detector section over the default section:

```julia
det = :DET_ID   # channel name
dsp_config_pd_det = merge(dsp_config_pd.default, get(dsp_config_pd, det, PropDict()))
config = DSPConfig(dsp_config_pd_det)
```

Keys of the detector section replace the keys of the default section, and nested sections are merged key by key. The layout expected by `DSPConfig(pd)` and the meaning of every field are described in [Waveforms and Configuration](@ref).

Which configuration type a function takes depends on the detector type:

| Detector type | Source of the configuration | Type passed to LegendDSP.jl | Functions |
|:--------------|:----------------------------|:----------------------------|:----------|
| HPGe | `dataprod_config(l200).dsp(filekey)` | [`DSPConfig`](@ref) | `dsp_icpc_compressed`, `dsp_decay_times` and all scans |
| Pulser and auxiliary channels | `dataprod_config(l200).dsp(filekey)`, listed under `additional_detectors` | [`DSPConfig`](@ref) | `dsp_puls_compressed`, `dsp_aux_compressed` |
| SiPM | `dataprod_config(l200).sipm(filekey)` | `PropDict` | `dsp_sipm_compressed`, `dsp_sg_sipm_optimization_compressed`, `dsp_sg_sipm_thresholds_compressed` |
| PMT | `dataprod_config(l200).pmt(filekey)` | `PropDict` | [`dsp_pmts`](@ref) |

For SiPM and PMT channels the merged `PropDict` is passed unchanged; the functions read the keys they need from it (see [SiPM and PMT Processing](@ref)).

The DSP section of the metadata also has subtrees that only the processors read and that are not part of `DSPConfig`: the `flt_optimization`, `aoe_optimization` and `pz` sections with the settings of the optimization fits and the quality selections, the lists `required_fltopt` and `required_aoeopt` with the filters whose parameters must exist before the DSP step runs, the selection expressions for the channels, and the mapping `additional_detectors`. The SiPM `optimization` section is passed to the SiPM scans as an argument.

## From Raw Data to DSP Output

The processors for calibration and physics data run these steps for each HPGe channel:

```julia
config = DSPConfig(merge(dsp_config_pd.default, get(dsp_config_pd, det, PropDict())))

raw_data = read_ldata(l200, DataTier(:raw), filekey, det)

dsp_data = dsp_icpc_compressed(raw_data, config, pars_tau[det].τ, pars_filter[det];
                               f_evaluate_qc = f_evaluate_qc)
```

`pars_tau[det].τ` is the decay time of the detector (see [Decay Time](@ref)). `pars_filter[det]` is the merge of the optimal energy filter parameters (`fltopt`) and current filter parameters (`aoeopt`) of the channel; see [Filter Optimization](@ref). Juleana.jl passes `missing` as the quality-cut evaluator, so `qc_label` is `-1`.

## Which Processor Calls Which Function

| Function | Processors | Call | Purpose |
|:---------|:-----------|:-----|:--------|
| `DSPConfig` | all HPGe, pulser and auxiliary processors, decay time, filter and A/E optimization | `DSPConfig(merge(default, detector section))` | Typed configuration of one channel |
| `dsp_icpc_compressed` | DSP for calibration and physics data | `(raw, config, τ, pars_filter; f_evaluate_qc)` | Main HPGe DSP |
| `dsp_decay_times` | decay time determination | `(wvfs, config)` | Decay time of each waveform of a gamma line; a truncated Gaussian is fitted to the distribution |
| `dsp_qc_flt_optimization_compressed` | decay time determination, filter optimization | `(wvfs, config, 400.0u"µs", f_evaluate_qc)` for the decay time, `(wvfs_pre, config, τ, f_evaluate_qc)` for the filters | Quick DSP for the waveform selection |
| `dsp_qdrift_flt_optimization` | filter optimization | `(wvfs_wdw, blmean_wdw, config, τ)` | Q-drift values for the flat-top fit |
| `dsp_{trap,cusp,zac}_rt_optimization` | filter optimization | `(wvfs_pre, config, τ; ft = ft_fixed)` | Noise grid over the rise times |
| `dsp_{trap,cusp,zac}_ft_optimization` | filter optimization | `(wvfs_pre, config, τ, rt)` | Energy grid over the flat-top times at the optimal rise time |
| `dsp_{sg,mwa}_optimization_compressed` | A/E optimization | `(wvfs_wdw, wvfs_pre, config, τ, pars_filter; f_evaluate_qc, presum_rate)` on waveforms of two calibration lines | A/E grid over the window lengths |
| `dsp_puls_compressed`, `dsp_aux_compressed` | DSP for physics data and for auxiliary channels | `(raw, config)` | Pulser, baseline and veto channels |
| [`dsp_pmts`](@ref) | DSP for physics data | `(raw, pmt_config)` | PMT channels |
| `dsp_sipm_compressed` | DSP for physics data | `(raw, sipm_config, pars_sipm[det])` | SiPM channels |
| `dsp_sg_sipm_optimization_compressed` | SiPM optimization | `(10000, wvfs, sipm_config, optimization_config)` | Trigger maxima and thresholds over the window lengths |
| `dsp_sg_sipm_thresholds_compressed` | SiPM optimization | `(wvfs, wl, sipm_config)` | Baseline and threshold distributions at the optimal window length |

The functions that are not in the table (`dsp_icpc`, `dsp_puls`, `dsp_aux`, `dsp_sipm`, `dsp_qc_flt_optimization`, `dsp_sg_optimization`, `dsp_mwa_optimization`, the primitives of [Waveform Parameters](@ref), the filters of [Filters](@ref) and the classifiers of [Machine-Learning Quality Cuts](@ref)) are not called by name from the processors. The processors work on the compressed raw data format, so they use the `_compressed` variants.

## Name-Based Dispatch

The processors do not hard-code the names of the optimization functions or of the auxiliary channel functions. They build a function name from a key of the configuration and look it up in the LegendDSP module:

```julia
function legend_dsp_function(name, config_key::AbstractString)
    sym = Symbol(name)
    isdefined(LegendDSP, sym) || throw(ArgumentError(
        "dataprod config key `$config_key` selects `$sym`, which LegendDSP does not define"))
    return getglobal(LegendDSP, sym)
end
```

An unknown name raises an `ArgumentError` that names the configuration key. The following table lists where the names come from and which functions of LegendDSP.jl each configuration can select.

| Configuration key | Name built by the processor | Valid values | Selected functions |
|:------------------|:----------------------------|:-------------|:-------------------|
| `dsp.additional_detectors.<channel>` | the string value, used as is | `"dsp_puls_compressed"`, `"dsp_aux_compressed"` | `dsp_puls_compressed`, `dsp_aux_compressed` (the uncompressed `dsp_puls` and `dsp_aux` have the same call signature) |
| `dsp.flt_optimization.e_filter.<filter>` | `dsp_<filter>_rt_optimization` and `dsp_<filter>_ft_optimization` | `trap`, `cusp`, `zac` | the six rise-time and flat-top scans |
| `dsp.aoe_optimization.aoe_filter.<filter>` | `dsp_<filter>_optimization_compressed` | `sg`, `mwa` | `dsp_sg_optimization_compressed`, `dsp_mwa_optimization_compressed` |
| `sipm.optimization.e_filter.<filter>` | `dsp_<filter>_sipm_optimization_compressed` and `dsp_<filter>_sipm_thresholds_compressed` | `sg` | `dsp_sg_sipm_optimization_compressed`, `dsp_sg_sipm_thresholds_compressed` |

The dispatch fixes more than the function name:

* For the energy filters the processor also reads the grid fields `e_grid_rt_<filter>` and `e_grid_ft_<filter>` of the `DSPConfig`, and it expects the scans to return the grid layout described in [Filter Optimization](@ref). The default configurations select `trap` and `cusp`; `zac` works when the configuration lists it and has the `e_grid_zac` entries.
* The keys `<filter>` of `aoe_filter` also name the entries of the optimal-parameter set `pars_filter` (`sg` and `mwa`, each with `wl`) that `dsp_icpc_compressed` reads.
* A new scan can be selected from the configuration by defining a function with a matching name in LegendDSP.jl and adding the settings it reads to the configuration.

The same pattern is used by the legacy dataflow, which looks the names up with `getfield(LegendDSP, Symbol(name))`.

## Typical Call Sequences

The determination of the parameters that the DSP step needs follows a fixed order. Each step stores its result as a parameter set that later steps read.

```julia
# 1. Decay time (pole-zero), per detector
dsp_qc  = dsp_qc_flt_optimization_compressed(wvfs, config, 400.0u"µs", f_evaluate_qc)
wvfs    = wvfs[findall(ljl_propfunc(qc_string).(dsp_qc))]     # quality selection
taus    = dsp_decay_times(wvfs, config)                        # then fit the distribution

# 2. Energy filter parameters, per detector and filter
dsp_qc  = dsp_qc_flt_optimization_compressed(wvfs_pre, config, τ, f_evaluate_qc)
qdrift  = dsp_qdrift_flt_optimization(wvfs_wdw, dsp_qc.blmean ./ presum_rate, config, τ)
for filter_type in keys(optimization_config.e_filter)
    rt_scan = legend_dsp_function("dsp_$(filter_type)_rt_optimization", "...")
    ft_scan = legend_dsp_function("dsp_$(filter_type)_ft_optimization", "...")
    enc_grid = rt_scan(wvfs_pre, config, τ; ft = optimization_config.e_filter[filter_type].ft_fixed)
    # fit_enc_sigmas(enc_grid, ...) -> optimal rt
    e_grid   = ft_scan(wvfs_pre, config, τ, rt_optimal)
    # fit_fwhm_ft(e_grid, ..., qdrift, ...) -> optimal ft
end

# 3. Current filter parameters, per detector and filter
for filter_type in keys(aoe_optimization_config.aoe_filter)
    f = legend_dsp_function("dsp_$(filter_type)_optimization_compressed", "...")
    dsp_dep = f(dep_wdw, dep_pre, config, τ, pars_filter; f_evaluate_qc, presum_rate)
    dsp_sep = f(sep_wdw, sep_pre, config, τ, pars_filter; f_evaluate_qc, presum_rate)
    # fit_sf_wl(dep.energy, dep.aoe, sep.energy, sep.aoe, grid) -> optimal wl
end

# 4. DSP of the calibration and physics data with all parameters
dsp_icpc_compressed(raw, config, τ, pars_filter; f_evaluate_qc)
```

The fits (`fit_enc_sigmas`, `fit_fwhm_ft`, `fit_sf_wl`) are part of LegendSpecFits.jl. The SiPM parameters are determined in a similar order with the SiPM scans: `dsp_sg_sipm_optimization_compressed` gives trigger maxima for each window length, `fit_sipm_wl` chooses the window length, `dsp_sg_sipm_thresholds_compressed` gives the threshold distributions at that window length, and the result is passed to `dsp_sipm_compressed` as `pars_sipm`.
