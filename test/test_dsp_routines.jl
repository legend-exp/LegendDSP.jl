using Test
using LegendDSP
using RadiationDetectorSignals
using RadiationDetectorDSP
using Unitful
using IntervalSets

function make_rise_waveforms()
    t = (0:8191) .* 16.0u"ns"
    signal = fill(1000.0, length(t))
    signal[3001:3125] .= range(1000.0, 11000.0; length=125)
    signal[3126:end] .= 11000.0 .* exp.(-((1:length(t)-3125) ./ (500e3 / 16))) .+ 1000.0
    ArrayOfRDWaveforms([RDWaveform(t, signal)])
end

function make_pulse_waveforms()
    t = (0:255) .* 16.0u"ns"
    signal = zeros(Float64, length(t))
    signal[101:110] .= 10.0
    ArrayOfRDWaveforms([RDWaveform(t, signal)])
end

@testset "dsp_routines" begin
    wvfs_rise = make_rise_waveforms()
    t0 = get_t0(wvfs_rise, 4.0; mintot=1500u"ns")[1]
    t50 = get_threshold(wvfs_rise, [0.5 * maximum(first(wvfs_rise).signal)]; mintot=1000u"ns")[1]
    @test 45u"µs" < t0 < 51u"µs"
    @test t0 < t50

    wvfs_pulse = make_pulse_waveforms()
    bl = 0u"ns" .. 80u"ns"
    trigs = get_triggers(wvfs_pulse, 5.0, bl; mintot=32u"ns")
    trigs_rev = get_triggers_reversed(wvfs_pulse, 5.0, bl; mintot=32u"ns")
    @test trigs.n == [1]
    @test trigs_rev.n == [1]
    @test trigs.intersect[1] < trigs_rev.intersect[1]

    flat = ArrayOfRDWaveforms([RDWaveform((0:127) .* 16.0u"ns", fill(1.0, 128))])
    qdrift = get_qdrift(wvfs_rise, [45.0u"µs"], 1.0u"µs":0.1u"µs":2.0u"µs")
    @test length(qdrift) == 1
    @test isfinite(ustrip(qdrift[1]))
end