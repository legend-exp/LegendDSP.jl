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

    # Edge case: Maximum at last sample of search window (ind_max == length(idxs))
    # Use short maxtot so the search window ends at the waveform end
    intflt_short = IntersectMaximum(mintot = 2Δt, maxtot = 5Δt)
    signal_max_at_last = zeros(Float64, n_samples)
    signal_max_at_last[end-4] = 0.3
    signal_max_at_last[end-3] = 0.5   # Crosses threshold here
    signal_max_at_last[end-2] = 0.6
    signal_max_at_last[end-1] = 0.8
    signal_max_at_last[end] = 1.0     # Maximum at very last sample
    wvf_max_at_last = RDWaveform(times, signal_max_at_last)

    res_max_last = intflt_short(wvf_max_at_last, 0.4)
    @test res_max_last.multiplicity == 1
    @test res_max_last.max[1] == 1.0  # No interpolation, just the value at the edge

    # Edge case: Intersection at last sample (threshold crossed near the end)
    signal_last_intersect = zeros(Float64, n_samples)
    signal_last_intersect[end-2] = 0.3
    signal_last_intersect[end-1] = 0.5   # Crosses threshold
    signal_last_intersect[end] = 0.6     # Still above, but need mintot samples
    wvf_last_intersect = RDWaveform(times, signal_last_intersect)

    # With mintot = 2Δt, we need 2 samples above threshold
    res_last_intersect = intflt(wvf_last_intersect, 0.4)
    @test res_last_intersect.multiplicity == 1
    @test res_last_intersect.x[1] > times[end-3]
end
