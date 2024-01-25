# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

import Test
import Aqua
import LegendDSP

Test.@testset "Aqua tests" begin
    Aqua.test_all(
        LegendDSP,
        ambiguities = false        
    )
end # testset
