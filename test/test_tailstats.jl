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

    # Noisy tail on a Float32 time axis far from zero, against a two-pass BigFloat fit of log(y)
    ts = range(60000f0u"ns", step = 16f0u"ns", length = 4000)
    t_ns = ustrip.(u"ns", collect(ts))
    ys = 2e4 .* exp.(-t_ns ./ 4e5) .+ 170 .* sin.(1.7 .* (1:length(ts)) .^ 2)
    tb, lb = big.(t_ns), log.(big.(ys))
    mt, ml = sum(tb) / length(tb), sum(lb) / length(lb)
    ref_sigma = sqrt(sum((lb .- ml) .^ 2) / length(lb))
    ref_τ = -sum((tb .- mt) .^ 2) / sum((tb .- mt) .* (lb .- ml))
    for (T, rtol) in ((Float64, 1e-10), (Float32, 1e-4))
        stats = tailstats(RDWaveform(ts, T.(ys)), first(ts), last(ts))
        @test stats.mean ≈ ml rtol = rtol
        @test stats.sigma ≈ ref_sigma rtol = rtol
        @test ustrip(u"ns", stats.τ) ≈ ref_τ rtol = rtol
    end
end
