using Test
using LegendDSP
using RadiationDetectorDSP
using RadiationDetectorSignals
using Unitful

@testset "Test WeightedSavitzkyGolayFilter" begin
    n = 600
    reltol = 1e-6
    t = range(0u"μs", 20u"μs", 2*n)
    signal = vcat(zeros(n), 10*ones(n))
    wf = RDWaveform(t, signal)

    # define filter parameters and filter
    flt = WeightedSavitzkyGolayFilter(length=1u"μs", degree=3, weightType=2)

    # apply filter to signal
    @test_nowarn flt(wf)
    @test isapprox(wf.signal[end], flt(wf).signal[end]; rtol=reltol)    
end

@testset "Test ModifiedSincFilter" begin
    n = 600
    reltol = 1e-6
    t = range(0u"μs", 20u"μs", 2*n)
    signal = vcat(zeros(n), 10*ones(n))
    wf = RDWaveform(t, signal)
    
    # define filter parameters and filter
    flt = ModifiedSincFilter(d=4, m=1u"μs")
    
    # apply filter to signal
    @test_nowarn flt(wf)
    @test isapprox(wf.signal[end], flt(wf).signal[end]; rtol=reltol)
end

@testset "Test WhittakerHendersonFilter" begin
    n = 600
    reltol = 1e-6
    t = range(0u"μs", 20u"μs", 2*n)
    signal = vcat(zeros(n), 10*ones(n))
    wf = RDWaveform(t, signal)

    # define filter parameters and filter
    flt = WhittakerHendersonFilter(p=3, λ=4)

    # apply filter to signal
    @test_nowarn flt(wf)
    @test isapprox(wf.signal[end], flt(wf).signal[end]; rtol=reltol)
end