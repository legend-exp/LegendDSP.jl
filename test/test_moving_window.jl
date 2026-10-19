using Test
using LegendDSP
using RadiationDetectorSignals
using RadiationDetectorDSP
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

@testset "MovingWindowMultiFilter against references" begin
    # moving mean over the last l samples, with samples before the start set to the first sample
    padded_mean(x, l) = [sum(x[max(m, 1)] for m in k-l+1:k) / l for k in eachindex(x)]
    times = (0:199) .* 16.0u"ns"
    mw(s, l) = MovingWindowFilter(l * step(times))(RDWaveform(times, s)).signal
    for (T, rtol) in ((Float64, 1e-10), (Float32, 1e-5)), l in (1, 2, 7, 36, 200, 250),
            x in (T.(randn(200)), T.(1000 .* (eachindex(times) .> 80) .+ randn(200)))
        out = MovingWindowMultiFilter(l * step(times))(RDWaveform(times, x)).signal
        @test eltype(out) == T
        @test out ≈ padded_mean(reverse(padded_mean(reverse(padded_mean(Float64.(x), l)), l)), l) rtol = rtol
        # same operations as MovingWindowFilter forward, backward, forward (defined for l ≤ length)
        if l <= length(x)
            @test out == mw(reverse(mw(reverse(mw(x, l)), l)), l)
        end
    end

    # integer input, and views as input and output
    xi = rand(-100:100, 200)
    @test MovingWindowMultiFilter(7 * step(times))(RDWaveform(times, xi)).signal ==
          MovingWindowMultiFilter(7 * step(times))(RDWaveform(times, float.(xi))).signal
    fi = RadiationDetectorDSP.fltinstance(MovingWindowMultiFilter(7 * step(times)), RadiationDetectorDSP.smplinfo(RDWaveform(times, zeros(200))))
    x, ybuf = randn(220), zeros(240)
    y = view(ybuf, 21:220)
    RadiationDetectorDSP.rdfilt!(y, fi, view(x, 11:210))
    @test y == RadiationDetectorDSP.rdfilt!(similar(x, 200), fi, x[11:210])
    @test all(iszero, ybuf[1:20]) && all(iszero, ybuf[221:240])
    @test_throws DimensionMismatch RadiationDetectorDSP.rdfilt!(zeros(199), fi, x[11:210])
end