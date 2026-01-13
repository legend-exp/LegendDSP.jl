using Test
using LegendDSP
using RadiationDetectorSignals
using Unitful

@testset "IntersectMaximum edge cases" begin
    Δt = 16u"ns"
    n_samples = 6200
    times = (0:n_samples-1) .* Δt
    intflt = IntersectMaximum(mintot = 2Δt, maxtot = 100Δt)

    # Intersection close to the beginning of the waveform
    signal_start = zeros(Float64, n_samples)
    signal_start[1] = 0.0
    signal_start[2] = 0.5
    signal_start[3] = 0.6
    signal_start[4] = 0.2
    wvf_start = RDWaveform(times, signal_start)

    res_start = intflt(wvf_start, 0.4)
    @test res_start.multiplicity >= 0
    @test length(res_start.x) == length(res_start.max)

    # Intersection close to the end of the waveform
    signal_end = zeros(Float64, n_samples)
    signal_end[end-3] = 0.0
    signal_end[end-2] = 0.5
    signal_end[end-1] = 0.6
    signal_end[end] = 0.2
    wvf_end = RDWaveform(times, signal_end)

    res_end = intflt(wvf_end, 0.4)
    @test res_end.multiplicity >= 0
    @test length(res_end.x) == length(res_end.max)
end
