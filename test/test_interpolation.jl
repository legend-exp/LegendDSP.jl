using Test
using LegendDSP
using RadiationDetectorSignals
using StructArrays: StructArray
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

    # A time axis with Float32 elements and a Float64 step, with one maximum at the window
    # edge and one inside the window
    times_mixed = StepRangeLen{typeof(0.0f0u"ns")}(0.0u"ns", 16.0u"ns", n_samples)
    wvf_edge = RDWaveform(times_mixed, Float32.(signal_start))
    wvf_inner = RDWaveform(times_mixed, Float32.(signal_mid))
    @inferred get_wvf_maximum(wvf_edge, times_mixed[1], times_mixed[54])
    @inferred get_wvf_maximum(wvf_inner, times_mixed[1], times_mixed[54])
    max_mixed = get_wvf_maximum.(ArrayOfRDWaveforms([wvf_edge, wvf_inner]), times_mixed[1], times_mixed[54])
    @test max_mixed isa StructArray
    @test max_mixed.t[1] == times_mixed[1]
    @test max_mixed.t[2] ≈ times_mixed[51]
    @test max_mixed.max[1] == 1.0f0
end

@testset "get_wvf_maximum index selection" begin
    times = (0:99) .* 16.0u"ns"

    # Ties: the first maximum in the window is used
    signal_tie = zeros(100)
    signal_tie[[19, 20, 21, 29, 30, 31]] .= [0.5, 1.0, 0.25, 0.25, 1.0, 0.5]
    m_tie = get_wvf_maximum(RDWaveform(times, signal_tie), times[10], times[40])
    @test m_tie.max == LegendDSP.extrema3points(0.5, 1.0, 0.25)
    @test times[19] < m_tie.t < times[21]

    # A NaN inside the window propagates to the result
    signal_nan = copy(signal_tie)
    signal_nan[25] = NaN
    @test isnan(get_wvf_maximum(RDWaveform(times, signal_nan), times[10], times[40]).max)
    @test !isnan(get_wvf_maximum(RDWaveform(times, signal_nan), times[26], times[40]).max)

    # Agreement with an argmax-based search on random signals and random windows
    function reference_wvf_maximum(X, Y, idxs)
        maxidx = idxs[argmax(Y[idxs])]
        if first(idxs) < maxidx < last(idxs)
            y1, y2, y3 = Y[maxidx-1], Y[maxidx], Y[maxidx+1]
            δ = (y1 - y3) / (2 * (y1 - 2y2 + y3))
            (max = LegendDSP.extrema3points(y1, y2, y3), t = X[maxidx] + step(X) * δ)
        else
            (max = Y[maxidx], t = X[maxidx] + step(X) * 0.0)
        end
    end
    for _ in 1:200
        Y = round.(randn(100); digits = 1)   # rounding produces ties
        i, j = minmax(rand(1:100), rand(1:100))
        @test get_wvf_maximum(RDWaveform(times, Y), times[i], times[j]) == reference_wvf_maximum(times, Y, i:j)
    end

    # A view of a longer signal gives the same result as a copy
    long_signal = randn(120)
    wvf_view = RDWaveform(times, view(long_signal, 11:110))
    wvf_copy = RDWaveform(times, long_signal[11:110])
    @test get_wvf_maximum(wvf_view, times[5], times[95]) == get_wvf_maximum(wvf_copy, times[5], times[95])
end
