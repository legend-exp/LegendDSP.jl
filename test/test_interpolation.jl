using Test
using LegendDSP
using RadiationDetectorSignals
using Unitful

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
    @test max_start.max >= 1.0   # maximum is at signal[1] = 1.0
    @test max_start.max < 1.1
    @test max_start.t == 0u"ns"

    # Maximum at the end of the search window
    signal_end = zeros(Float64, n_samples)
    signal_end[end-3] = 0.2
    signal_end[end-2] = 0.5
    signal_end[end-1] = 0.8
    signal_end[end] = 1.0
    wvf_end = RDWaveform(times, signal_end)

    max_end = get_wvf_maximum(wvf_end, times[end-4], times[end])
    @test max_end.max >= 1.0   # maximum is at signal[end] = 1.0
    @test max_end.max < 1.1
    @test max_end.t == times[end]

    # Maximum in the middle (uses parabola interpolation)
    signal_mid = zeros(Float64, n_samples)
    signal_mid[50] = 0.5
    signal_mid[51] = 1.0
    signal_mid[52] = 0.5
    wvf_mid = RDWaveform(times, signal_mid)

    max_mid = get_wvf_maximum(wvf_mid, times[48], times[54])
    @test max_mid.max >= 1.0   # interpolated maximum should be at least 1.0
    @test max_mid.max < 1.2    # but not unreasonably high

    # Quadratic interpolation of a maximum between two samples
    signal_time = @. -(ustrip(u"ns", times) - 805)^2
    wvf_time = RDWaveform(times, signal_time)
    @test get_wvf_maximum(wvf_time, times[48], times[54]).t ≈ 805u"ns"

    # Broadcasting over waveform StructArrays returns a dense vector of times
    times_float = (0:n_samples-1) .* 16.0u"ns"
    wvf_time_float = RDWaveform(times_float, signal_time)
    wvfs_time = ArrayOfRDWaveforms([wvf_time_float, wvf_time_float])
    max_times = get_wvf_maximum.(wvfs_time, times_float[48], times_float[54]).t
    @test max_times isa Vector
    @test max_times ≈ fill(805u"ns", 2)
end
