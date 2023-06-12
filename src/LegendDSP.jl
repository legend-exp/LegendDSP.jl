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

import Adapt


include("tailstats.jl")

end # module
