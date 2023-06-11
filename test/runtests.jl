# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

import Test

Test.@testset "Package LegendDSP" begin
    include("test_aqua.jl")
    # include("test_SOMETHING.jl")
    include("test_docs.jl")
    isempty(Test.detect_ambiguities(LegendDSP))
end # testset
