using Test
using LegendDSP
using RadiationDetectorDSP
using RadiationDetectorSignals
using Unitful

@testset "Test WhittakerHendersonFilter" begin
    n = 600
    noise = 1.
    t = range(0u"μs", 20u"μs", 2*n)
    signal = vcat(zeros(n), 10*ones(n)) + (noise*rand(2*n) .- noise/2)
    wf = RDWaveform(t, signal)

    # define filter parameters and filter
    flt = WhittakerHendersonFilter(p=3, λ=4)

    # apply filter to signal
    @test_nowarn wf_new = flt(wf)
end