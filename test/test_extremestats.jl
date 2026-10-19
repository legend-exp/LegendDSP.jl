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

    @testset "Test index selection" begin
        times = (0:99) .* 16.0u"ns"

        # Ties: the first minimum and the first maximum are used
        signal_tie = zeros(100)
        signal_tie[[20, 30]] .= 1.0
        signal_tie[[25, 35]] .= -1.0
        es_tie = extremestats(RDWaveform(times, signal_tie))
        @test es_tie == (min = -1.0, max = 1.0, tmin = times[25], tmax = times[20])

        # A constant signal has both extrema at the first sample of the window
        es_const = extremestats(RDWaveform(times, fill(2.0, 100)), times[10], times[90])
        @test es_const == (min = 2.0, max = 2.0, tmin = times[10], tmax = times[10])

        # A NaN in the window is returned as minimum and maximum, at its position
        signal_nan = copy(signal_tie)
        signal_nan[40] = NaN
        es_nan = extremestats(RDWaveform(times, signal_nan))
        @test isnan(es_nan.min) && isnan(es_nan.max)
        @test es_nan.tmin == es_nan.tmax == times[40]

        # Agreement with findmin/findmax on random signals (with ties and signed zeros) and windows
        for _ in 1:200
            Y = round.(randn(100); digits = 1) .* rand([-1.0, 1.0], 100)
            i, j = minmax(rand(1:100), rand(1:100))
            (minv, imin), (maxv, imax) = findmin(Y[i:j]), findmax(Y[i:j])
            @test extremestats(RDWaveform(times, Y), times[i], times[j]) ===
                (min = minv, max = maxv, tmin = times[i + imin - 1], tmax = times[i + imax - 1])
        end

        # A view of a longer signal gives the same result as a copy
        long_signal = randn(120)
        @test extremestats(RDWaveform(times, view(long_signal, 11:110)), times[5], times[95]) ==
            extremestats(RDWaveform(times, long_signal[11:110]), times[5], times[95])
    end
end
