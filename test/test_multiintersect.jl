using Test
using LegendDSP
using RadiationDetectorSignals
using RadiationDetectorDSP
using Unitful

@testset "Test MultiIntersect filter" begin
    @testset "MultiIntersect similar to Intersect" begin
        wvf = RDWaveform((1:100)u"s", collect(1:100))

        flt_mult = MultiIntersect(threshold_ratios=[0.5], mintot=1u"s")
        flt = Intersect(mintot=1u"s")
        @test isapprox(flt(wvf, 0.5*100).x, flt_mult(wvf)[1])
    end
    @testset "Linear Waveform" begin
        
        wvf = RDWaveform((1:100)u"s", collect(1:100))

        @testset "MultiIntersect(1u\"s\")" begin
            ratios = 0.1:0.1:0.9 |> collect
            flt = MultiIntersect(threshold_ratios=ratios, mintot=1u"s")

            res = flt(wvf)
            correct_res = collect(10:10:90)*u"s"
            @test isapprox(res, correct_res)
        end
    end
end