# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    LegendDSP

Template for Julia packages.
"""
module LegendDSP

using LinearAlgebra
using Random

using ArgCheck
using ArraysOfArrays
using DocStringExtensions
using FillArrays
using FunctionChains
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
include("types.jl")
include("utils.jl")
include("derivative.jl")
include("dsp_decaytime.jl")
include("dsp_filter_optimization.jl")
include("dsp_aoe_optimization.jl")
include("dsp_icpc.jl")
include("saturation.jl")
include("interpolation.jl")
include("intersect_maximum.jl")
include("dsp_sipm.jl")
include("dsp_sipm_optimization.jl")
include("dsp_routines.jl")
include("dsp_puls.jl")
include("haar_filter.jl")
include("ml.jl")
include("dsp_ml_routines.jl")
include("multi_intersect.jl")
include("moving_window_multi.jl")
include("alternative_filters.jl")

end # module
