# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

using Test
using LegendDSP
using RadiationDetectorSignals
using Unitful

@testset "Test extremestats" begin
    
    @testset "Test single waveform" begin
        wf = RDWaveform((0:1:360)u"ns", sind.(0:360))

        es_full = extremestats(wf)
        @test es_full.min  == -1.0
        @test es_full.max  ==  1.0
        @test es_full.tmin ==  270u"ns"
        @test es_full.tmax ==   90u"ns"

        es_part = extremestats(wf, 0u"ns", 180u"ns")
        @test es_part.min  == 0.0
        @test es_part.max  == 1.0
        @test es_part.tmin ==   0u"ns"
        @test es_part.tmax ==  90u"ns"

        es_part2 = extremestats(wf, 135u"ns", 225u"ns")
        @test es_part2.min  == -sqrt(0.5)
        @test es_part2.max  ==  sqrt(0.5)
        @test es_part2.tmin == 225u"ns"
        @test es_part2.tmax == 135u"ns"
    end
    
    @testset "Test single array" begin
        wf = sind.(1:360)
        
        es_full = extremestats(wf)
        @test es_full.min  == -1.0
        @test es_full.max  ==  1.0
        @test es_full.tmin ==  270
        @test es_full.tmax ==   90

        es_part = extremestats(wf, 180, 360)
        @test es_part.min  == -1.0
        @test es_part.max  ==  0.0
        @test es_part.tmin ==  270
        @test es_part.tmax ==  180

        es_part2 = extremestats(wf, 135, 225)
        @test es_part2.min  == -sqrt(0.5)
        @test es_part2.max  ==  sqrt(0.5)
        @test es_part2.tmin == 225
        @test es_part2.tmax == 135
    end
end
