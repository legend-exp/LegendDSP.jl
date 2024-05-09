using Test
using LegendDSP
using RadiationDetectorSignals
using Unitful

@testset "Test MultiIntersect filter" begin
    @testset "Linear Waveform" begin
        
        wvf = RDWaveform((1:100)u"s", collect(1:100))

        @testset "MultiIntersect(1u\"s\")" begin
            flt = MultiIntersect(1u"s")

            res = flt(wvf)
            correct_res = collect(10:10:90)*u"s"
            @test isapprox(res, correct_res)
        end
    end
end