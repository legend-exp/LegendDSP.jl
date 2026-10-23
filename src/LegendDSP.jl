# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    LegendDSP

Digital signal processing for LEGEND waveforms.

The package provides

- DSP chains that reduce tables of waveforms to tables of waveform parameters: [`dsp_icpc`](@ref) for HPGe detectors, [`dsp_puls`](@ref) for pulser channels, [`dsp_aux`](@ref) for auxiliary channels, [`dsp_sipm`](@ref) for SiPMs and [`dsp_pmts`](@ref) for PMTs, each with a `_compressed` variant where compressed data exists;
- grid scans over filter parameters (`dsp_*_optimization`) that provide the input for choosing optimal filter settings;
- functions for single waveform parameters ([`get_t0`](@ref), [`get_threshold`](@ref), [`get_a`](@ref), [`get_qdrift`](@ref), [`get_lq`](@ref), [`get_triggers`](@ref)) and waveform statistics ([`tailstats`](@ref), [`extremestats`](@ref), [`thresholdstats`](@ref), [`saturation`](@ref));
- filters for `RadiationDetectorDSP.jl` waveforms, such as [`MovingWindowMultiFilter`](@ref), [`WeightedSavitzkyGolayFilter`](@ref) and [`HaarAveragingFilter`](@ref);
- support-vector-machine quality cuts ([`get_qc_ml_func`](@ref), [`get_qc_classifier`](@ref));
- the configuration type [`DSPConfig`](@ref).
"""
module LegendDSP

using LinearAlgebra
using Random
using Statistics

using ArgCheck
using ArraysOfArrays
using DocStringExtensions
using FillArrays
using IntervalSets
using RadiationDetectorDSP
using RadiationDetectorSignals
using Unitful
using PropDicts
using Tables
using TypedTables
using LIBSVM
using LegendDataTypes

import Adapt
import RadiationDetectorDSP: fltinstance, rdfilt!, flt_output_length, 
    flt_input_smpltype, flt_output_smpltype, _smpllen, smplinfo, SamplingInfo,
    AbstractRadSigFilterInstance, LinearFiltering, _filterlen, _floattype

using Unitful: RealOrRealQuantity as RealQuantity

include("tailstats.jl")
include("extremestats.jl")
include("thresholdstats.jl")
include("types.jl")
include("utils.jl")
include("derivative.jl")
include("dsp_decaytime.jl")
include("dsp_filter_optimization.jl")
include("dsp_icpc.jl")
include("saturation.jl")
include("interpolation.jl")
include("intersect_maximum.jl")
include("dsp_sipm.jl")
include("dsp_sipm_optimization.jl")
include("dsp_pmts.jl")
include("dsp_routines.jl")
include("dsp_aux.jl")
include("haar_filter.jl")
include("ml.jl")
include("dsp_ml_routines.jl")
include("multi_intersect.jl")
include("moving_window_multi.jl")
include("alternative_filters.jl")
include("timeaxis.jl")

end # module
