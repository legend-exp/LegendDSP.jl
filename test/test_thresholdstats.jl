using Test
using LegendDSP
using RadiationDetectorSignals
using Unitful
using Statistics

@testset "thresholdstats_mad" begin
    Δt = 16u"ns"

    # constant signal → MAD should be zero
    @testset "constant signal" begin
        t = (0:99) .* Δt
        sig = fill(5.0, 100)
        wvf = RDWaveform(t, sig)
        result = thresholdstats_mad(wvf, -10.0, 10.0)
        @test result ≈ 0.0 atol=1e-10
    end

    # symmetric signal: known MAD value
    @testset "symmetric signal" begin
        t = (0:99) .* Δt
        # median=0, all |deviations|=1, so MAD=1, scaled=1.4826
        sig = vcat(fill(-1.0, 50), fill(1.0, 50))
        wvf = RDWaveform(t, sig)
        result = thresholdstats_mad(wvf, -5.0, 5.0)
        @test result ≈ 1.4826 atol=1e-10
    end

    # robustness: outlier pulse should not dominate MAD
    @testset "outlier robustness" begin
        t = (0:999) .* Δt
        sig = zeros(1000)
        # inject large outlier pulse in a few samples
        sig[500:510] .= 1000.0
        wvf = RDWaveform(t, sig)
        result = thresholdstats_mad(wvf, -Inf, Inf)
        # MAD should remain small, dominated by the majority of zero samples
        @test result < 1.0
    end

    # empty filter: all samples outside bounds → return 0
    @testset "empty filter" begin
        t = (0:99) .* Δt
        sig = fill(5.0, 100)
        wvf = RDWaveform(t, sig)
        result = thresholdstats_mad(wvf, 10.0, 20.0)
        @test result ≈ 0.0 atol=1e-10
    end

    # unitful signal
    @testset "unitful signal" begin
        t = (0:99) .* Δt
        sig = fill(5.0u"mV", 100)
        wvf = RDWaveform(t, sig)
        result = thresholdstats_mad(wvf, -10.0u"mV", 10.0u"mV")
        @test result ≈ 0.0u"mV" atol=1e-10u"mV"
    end

    # plain array (non-waveform)
    @testset "plain array" begin
        sig = vcat(fill(-1.0, 50), fill(1.0, 50))
        result = thresholdstats_mad(sig, -5.0, 5.0)
        @test result ≈ 1.4826 atol=1e-10
    end
end
