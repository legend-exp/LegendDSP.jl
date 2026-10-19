# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

using Test
using LegendDSP
using RadiationDetectorSignals
using Unitful

@testset "tailstats" begin
    times = range(0.0f0u"ns", step = 16.0f0u"ns", length = 1000)
    for T in (Float32, Float64)
        wvf = RDWaveform(times, T.(1000 .* exp.(-ustrip.(u"ns", times) ./ 5000)))
        stats = @inferred tailstats(wvf, 1000u"ns", 10000u"ns")
        @test stats.τ ≈ 5000u"ns" rtol = 1e-3
        @test stats.sigma isa Float64

        # a window with non-positive samples gives zeros of the same types
        wvf_neg = RDWaveform(times, fill(-one(T), length(times)))
        stats_neg = @inferred tailstats(wvf_neg, 1000u"ns", 10000u"ns")
        @test typeof(stats_neg) == typeof(stats)
        @test iszero(stats_neg.mean) && iszero(stats_neg.sigma) && iszero(stats_neg.τ)
    end
end
