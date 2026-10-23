# LegendDSP.jl

LegendDSP.jl is a Julia package for the digital signal processing of detector waveforms of the [LEGEND](https://legend-exp.org) experiment. It contains the processing chains for high-purity germanium (HPGe) detectors, pulser and auxiliary channels, silicon photomultipliers (SiPM) and photomultiplier tubes (PMT), and the functions that scan filter parameters on sets of waveforms. It is part of the LEGEND Julia analysis software.

The processing chains take the table of raw waveforms of one channel and return a table with one row per event: baseline and tail statistics, timing points, energy estimators from trapezoidal, CUSP and ZAC filters, the maximum current used for the A/E pulse-shape parameter, charge-drift and late-charge parameters, pile-up and saturation counts, and an optional label from a machine-learning quality cut. The configuration of the HPGe chains is held in a typed [`DSPConfig`](@ref) that is built from the production metadata.

## Place in the LEGEND Julia Stack

LegendDSP.jl builds on two general-purpose packages and is used by the production dataflow:

* [RadiationDetectorSignals.jl](https://github.com/JuliaPhysics/RadiationDetectorSignals.jl) provides the waveform types (`RDWaveform`, `ArrayOfRDWaveforms`) with Unitful time axes.
* [RadiationDetectorDSP.jl](https://github.com/JuliaPhysics/RadiationDetectorDSP.jl) provides the generic filters and signal statistics (trapezoidal, CUSP, ZAC, Savitzky-Golay, CR and inverse CR filters, `signalstats`, `Intersect`, signal estimators). LegendDSP.jl adds filters and parameters specific to LEGEND data and combines them into chains.
* [Juleana.jl](https://github.com/legend-exp/Juleana.jl) is the production dataflow. Its processors read the raw data tier, call the functions of this package, and write the DSP tier. The energy calibration and pulse-shape fits that follow are in [LegendSpecFits.jl](https://github.com/legend-exp/LegendSpecFits.jl).

## Installation

LegendDSP.jl is registered in the [LEGEND Julia registry](https://github.com/legend-exp/LegendJuliaRegistry). After adding the registry once, the package can be installed with the Julia package manager:

```julia
using Pkg
pkg"registry add General https://github.com/legend-exp/LegendJuliaRegistry.git"
Pkg.add("LegendDSP")
```

To work on the package itself, use `Pkg.develop("LegendDSP")`.

## Overview

The manual pages describe the concepts and the workflow of each part of the package:

* [Waveforms and Configuration](@ref) describes the waveform types, the uncompressed and compressed raw data layouts, and the fields of [`DSPConfig`](@ref).
* [HPGe DSP Chains](@ref) lists the processing steps and all output columns of [`dsp_icpc`](@ref), [`dsp_icpc_compressed`](@ref) and the pulser and auxiliary chains.
* [Filters](@ref) covers the filter types of the package: the moving-window and Savitzky-Golay filters, the alternative smoothing filters, the Haar filter, the derivative filter and the time-axis filter.
* [Waveform Parameters](@ref) covers the functions that extract single parameters from waveforms (t0, threshold crossings, current, Q-drift, LQ and triggers) and the statistics on tails, extrema, thresholds and saturation.
* [Decay Time](@ref) explains how the decay time of the preamplifier signal is determined and how it enters the pole-zero correction.
* [Filter Optimization](@ref) covers the scans of rise time, flat-top time and current-filter window length.
* [SiPM and PMT Processing](@ref) covers the SiPM and PMT chains and the SiPM threshold and window-length scans.
* [Machine-Learning Quality Cuts](@ref) covers the support-vector-machine quality cut and its interface to the chains.
* [Use in Juleana.jl](@ref) describes how the production dataflow builds configurations and selects functions of this package by name.

The tutorials run complete examples on synthetic data: [Energy Reconstruction](@ref) and [A/E Pulse-Shape Parameter](@ref). The [API](@ref) page lists all documented functions and types.

## Related Packages

* [RadiationDetectorSignals.jl](https://github.com/JuliaPhysics/RadiationDetectorSignals.jl) and [RadiationDetectorDSP.jl](https://github.com/JuliaPhysics/RadiationDetectorDSP.jl) provide the waveform types and the generic filters on which this package builds.
* [LegendDataTypes.jl](https://github.com/legend-exp/LegendDataTypes.jl) provides the data types of the LEGEND file formats and `decode_data`, which decodes encoded waveform columns.
* [LegendDataManagement.jl](https://github.com/legend-exp/LegendDataManagement.jl) reads LEGEND data and metadata, including the DSP configuration.
* [LegendSpecFits.jl](https://github.com/legend-exp/LegendSpecFits.jl) fits the parameters produced by the chains, for example the energy calibration, the A/E cut, and the optimal filter parameters from the scans.
* [LegendMakie.jl](https://github.com/legend-exp/LegendMakie.jl) provides plot recipes for LEGEND data.
* [Juleana.jl](https://github.com/legend-exp/Juleana.jl) is the production dataflow that applies the functions of this package to LEGEND data.
