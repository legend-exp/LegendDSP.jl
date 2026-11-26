using Test
using LegendDSP
using RadiationDetectorSignals
using Unitful

@testset "IntersectExtrema" begin
    Δt = 16u"ns"
    n = 400
    times = (0:n-1) .* Δt
    bump(i0, σ) = [exp(-0.5 * ((i - i0) / σ)^2) for i in 1:n]
    intflt = IntersectExtrema(mintot = 3Δt)

    # two maxima and one minimum in separate regions
    y = 2.0 .* bump(100, 5.0) .+ 0.7 .* bump(220, 5.0) .- 1.2 .* bump(300, 5.0)
    wf = RDWaveform(times, y)
    res = intflt(wf, 0.3, 0.3)
    @test length(res.t) == length(res.amp) == length(res.t_start) == length(res.t_end) == 3
    @test isapprox(res.t, [99, 219, 299] .* Δt; atol = 0.1Δt)
    @test isapprox(res.amp, [2.0, 0.7, -1.2]; atol = 0.01)
    @test res.t_start[1] == times[findfirst(>(0.3), y)]
    @test res.t_end[1] == times[findlast(>(0.3), y[1:150]) + 1]
    @test res.t_start[3] == times[findfirst(<(-0.3), y)]
    @test all(res.t_start .< res.t .< res.t_end)
    @test eltype(res.t) == eltype(res.t_start) == float(eltype(wf.time))
    @test isempty(intflt(wf, 3.0, 0.3).t)

    # two maxima in one region: both kept if the valley between them is deeper than `valley`, else the larger one
    y2 = bump(100, 5.0) .+ 0.8 .* bump(118, 5.0)
    res2 = intflt(RDWaveform(times, y2), 0.3, 0.3)
    @test isapprox(res2.t, [99, 117] .* Δt; atol = 0.5Δt)
    @test isapprox(res2.amp, [1.0, 0.8]; atol = 0.01)
    @test res2.t_start[1] == res2.t_start[2]
    @test res2.t_end[1] == res2.t_end[2]
    res2_merged = intflt(RDWaveform(times, y2), 0.3, 0.6)
    @test length(res2_merged.t) == 1
    @test res2_merged.amp[1] ≈ 1.0 atol = 0.01
    @test res2_merged.t_start[1] == res2.t_start[1]

    # a flat step on a rising edge is no extremum, a flat top is one
    @test IntersectExtrema(mintot = 3)([0.0, 1.0, 2.0, 2.0, 3.0, 1.0, 0.0], 0.5, 0.0).amp ≈ [3 + 1 / 24]
    @test length(IntersectExtrema(mintot = 3)([0.0, 1.0, 3.0, 3.0, 3.0, 1.0, 0.0], 0.5, 0.0).t) == 1

    # a region shorter than mintot is not counted
    y_short = zeros(n)
    y_short[200:201] .= 1.0
    @test isempty(intflt(RDWaveform(times, y_short), 0.3, 0.3).t)
    @test length(IntersectExtrema(mintot = 2Δt)(RDWaveform(times, y_short), 0.3, 0.3).t) == 1

    # regions that touch the first or the last sample are not counted
    y_edge = zeros(n)
    y_edge[1:10] .= 1.0
    y_edge[end-9:end] .= 1.0
    @test isempty(intflt(RDWaveform(times, y_edge), 0.3, 0.3).t)

    @test_throws ArgumentError intflt(wf, -0.3, 0.3)
    @test_throws ArgumentError intflt(wf, 0.3, -0.3)
    @test_throws ArgumentError intflt(wf, NaN, 0.3)
end
