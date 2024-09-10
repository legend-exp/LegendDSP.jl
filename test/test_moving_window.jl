using Test
using LegendDSP
using RadiationDetectorSignals
using Unitful

@testset "Test MovingWindowFilter" begin
    T = Float64
    sig = vcat(zeros(T, 5), ones(T, 5))
    res = [zeros(T, 5)..., 0.5, 1., 1., 1., 1.]
    wvf = RDWaveform(range(0u"μs", step = 32u"μs", length = 10), sig)

    flt = MovingWindowFilter(64u"μs")
    wvf_out = flt(wvf)
    @test isapprox(wvf_out.signal, res)
end
@testset "Test MovingWindowMultiFilter" begin
    T = Float64
    sig = vcat(zeros(T, 5), ones(T, 5))
    res = [0., 0., 0., 0., 0.125, 0.5, 0.875, 1., 1., 1.]
    wvf = RDWaveform(range(0u"μs", step = 32u"μs", length = 10), sig)

    flt = MovingWindowMultiFilter(64u"μs")
    wvf_out = flt(wvf)
    @test isapprox(wvf_out.signal, res)
end