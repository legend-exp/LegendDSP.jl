using Test
using LegendDSP
using RadiationDetectorSignals
using Unitful
using Random

@testset "madstats" begin
    Δt = 16u"ns"
    t = (0:999) .* Δt

    # median 0, all deviations 1: sigma 1.4826
    sig = vcat(fill(-1.0, 500), fill(1.0, 500))
    res = madstats(RDWaveform(t, sig), -5.0, 5.0)
    @test res.median ≈ 0.0 atol = 1e-10
    @test res.sigma ≈ 1.4826 atol = 1e-10
    @test madstats(sig, -5.0, 5.0) == res

    # an offset and a pulse on top of noise: the median is the offset, sigma the noise width
    noise = 0.1 .* randn(MersenneTwister(1), 1000)
    sig_pulse = 5.0 .+ noise
    sig_pulse[500:520] .+= 100.0
    res_pulse = madstats(RDWaveform(t, sig_pulse))
    @test res_pulse.median ≈ 5.0 atol = 0.02
    @test res_pulse.sigma ≈ 0.1 atol = 0.02
    @test res_pulse.sigma == thresholdstats_mad(RDWaveform(t, sig_pulse))

    # only the samples between the bounds count
    sig_two = vcat(fill(0.0, 400), fill(10.0, 600))
    @test madstats(sig_two).median == 10.0
    @test madstats(sig_two, -1.0, 1.0).median == 0.0

    # no sample between the bounds
    res_empty = madstats(RDWaveform(t, fill(5.0, 1000)), 10.0, 20.0)
    @test res_empty.median == 0.0
    @test res_empty.sigma == 0.0

    # unitful signal
    res_unit = madstats(RDWaveform(t, fill(5.0u"mV", 1000)), -10.0u"mV", 10.0u"mV")
    @test res_unit.median == 5.0u"mV"
    @test res_unit.sigma == 0.0u"mV"
end

@testset "sigmaclipstats" begin
    # integer ADC values around 3750.3: the mean resolves the level below one count, the median does not
    rng = MersenneTwister(2)
    sig = round.(3750.3 .+ 1.3 .* randn(rng, 2000))
    res = sigmaclipstats(sig, 3.0)
    @test res.mean ≈ 3750.3 atol = 0.1
    @test res.sigma ≈ 1.4826 atol = 1e-10

    # samples farther than n_σ noise widths from the median do not enter the mean
    sig_pulse = copy(sig)
    sig_pulse[1000:1020] .+= 100.0
    @test sigmaclipstats(sig_pulse, 3.0).mean ≈ res.mean atol = 0.02
    @test sigmaclipstats(sig_pulse, 1000.0).mean > res.mean + 0.5

    # waveform and unitful signal
    t = (0:1999) .* 16u"ns"
    @test sigmaclipstats(RDWaveform(t, sig), 3.0) == res
    res_unit = sigmaclipstats(RDWaveform(t, sig .* 1.0u"mV"), 3.0)
    @test res_unit.mean ≈ res.mean * u"mV"
    @test res_unit.sigma ≈ res.sigma * u"mV"
end
