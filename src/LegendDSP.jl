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

import Adapt


include("tailstats.jl")
include("utils.jl")
include("decaytime.jl")
include("filter_optimization.jl")

@static if !isdefined(Base, :get_extension)
    using Requires
    include("../ext/LegendDSPRecipesBaseExt.jl")
end

end # module
