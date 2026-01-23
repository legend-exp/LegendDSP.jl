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
    @test res_start.multiplicity == 1
    @test length(res_start.x) == length(res_start.max)
    @test length(res_start.x) == 1
    @test res_start.x[1] > 0u"ns"
    @test res_start.x[1] < 48u"ns"  # intersection should be between sample 2 and 3
    @test res_start.max[1] >= 0.6   # maximum should be at least 0.6 (interpolated value can be higher)
    @test res_start.max[1] < 0.7    # but not unreasonably high

    # Intersection close to the end of the waveform
    signal_end = zeros(Float64, n_samples)
    signal_end[end-3] = 0.0
    signal_end[end-2] = 0.5
    signal_end[end-1] = 0.6
    signal_end[end] = 0.2
    wvf_end = RDWaveform(times, signal_end)

    res_end = intflt(wvf_end, 0.4)
    @test res_end.multiplicity == 1
    @test length(res_end.x) == length(res_end.max)
    @test length(res_end.x) == 1
    @test res_end.x[1] > times[end-3]
    @test res_end.x[1] < times[end]
    @test res_end.max[1] >= 0.6   # maximum should be at least 0.6 (interpolated value can be higher)
    @test res_end.max[1] < 0.7    # but not unreasonably high
end

@testset "get_wvf_maximum edge cases" begin
    Δt = 16u"ns"
    n_samples = 100
    times = (0:n_samples-1) .* Δt

    # Maximum at the beginning of the search window
    signal_start = zeros(Float64, n_samples)
    signal_start[1] = 1.0
    signal_start[2] = 0.8
    signal_start[3] = 0.5
    signal_start[4] = 0.2
    wvf_start = RDWaveform(times, signal_start)

    max_start = get_wvf_maximum(wvf_start, 0u"ns", 64u"ns")
    @test max_start >= 1.0   # maximum is at signal[1] = 1.0
    @test max_start < 1.1

    # Maximum at the end of the search window
    signal_end = zeros(Float64, n_samples)
    signal_end[end-3] = 0.2
    signal_end[end-2] = 0.5
    signal_end[end-1] = 0.8
    signal_end[end] = 1.0
    wvf_end = RDWaveform(times, signal_end)

    max_end = get_wvf_maximum(wvf_end, times[end-4], times[end])
    @test max_end >= 1.0   # maximum is at signal[end] = 1.0
    @test max_end < 1.1

    # Maximum in the middle (uses parabola interpolation)
    signal_mid = zeros(Float64, n_samples)
    signal_mid[50] = 0.5
    signal_mid[51] = 1.0
    signal_mid[52] = 0.5
    wvf_mid = RDWaveform(times, signal_mid)

    max_mid = get_wvf_maximum(wvf_mid, times[48], times[54])
    @test max_mid >= 1.0   # interpolated maximum should be at least 1.0
    @test max_mid < 1.2    # but not unreasonably high
end
