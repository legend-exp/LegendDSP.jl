using Test
using LegendDSP
using RadiationDetectorSignals
using Unitful
using Random

# a shaped response: Gaussian bump of unit peak, σ in samples
bump(n, i0, σ) = [exp(-0.5 * ((i - i0) / σ)^2) for i in 1:n]

@testset "CenteredConvolutionFilter" begin
    n = 200
    wf = RDWaveform(range(0u"ns", step = 16u"ns", length = n), bump(n, 80, 6.0))
    # a unit tap at lag 0 reproduces the signal; a tap at lag +2 delays it by two samples (y[i] = Σ h[lag] x[i − lag],
    # the convention of the kernel design in LegendSpecFits)
    @test CenteredConvolutionFilter([0.0, 0.0, 1.0, 0.0, 0.0])(wf).signal ≈ wf.signal
    out = CenteredConvolutionFilter([0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0])(wf)
    @test out.time == wf.time
    @test isapprox(out.signal[3:n], wf.signal[1:n-2]; atol = 1e-12)
    @test isapprox(out.signal[1:2], zeros(2); atol = 1e-12)
    # a two-sided kernel against the direct sum, zero beyond the ends
    h = [0.1, -0.3, 1.0, -0.3, 0.1]
    direct = [sum(h[l] * (1 <= i + 3 - l <= n ? wf.signal[i+3-l] : 0.0) for l in 1:5) for i in 1:n]
    @test CenteredConvolutionFilter(h)(wf).signal ≈ direct
    # a kernel longer than the signal (as for the template response in dsp_sipm)
    hl = sin.(1:401) ./ 50
    directl = [sum(hl[l] * (1 <= i + 201 - l <= n ? wf.signal[i+201-l] : 0.0) for l in 1:401) for i in 1:n]
    @test CenteredConvolutionFilter(hl)(wf).signal ≈ directl
    # one instance for many waveforms, as in a broadcast
    @test CenteredConvolutionFilter(h).([wf, wf])[2].signal ≈ direct
    @test_throws ArgumentError CenteredConvolutionFilter([1.0, 0.0])(wf)
end

@testset "extrema_over_threshold" begin
    n = 400
    rng = MersenneTwister(3)
    y = 2.0 .* bump(n, 100, 5.0) .+ 0.7 .* bump(n, 220, 5.0) .- 1.2 .* bump(n, 300, 5.0) .+ 0.002 .* randn(rng, n)
    wf = RDWaveform(range(0u"ns", step = 16u"ns", length = n), y)
    sites = extrema_over_threshold(wf, 0.3)
    @test length(sites.x) == 3
    @test isapprox(sites.x, [99, 219, 299] .* 16.0u"ns"; atol = 16.0u"ns")
    @test isapprox(sites.max, [2.0, 0.7, -1.2]; atol = 0.05)
    @test all(sites.x_tot .> 0u"ns")
    @test sites.x_tot[1] > sites.x_tot[2]        # the wider excursion of the larger bump
    # plain samples work too and nothing below threshold is a site
    @test length(extrema_over_threshold(y, 3.0).x) == 0
    @test length(extrema_over_threshold(y, 0.3).x) == 3
    # two photons whose excursions merge: both maxima are sites, sharing the time over threshold
    y2 = bump(n, 100, 5.0) .+ 0.8 .* bump(n, 115, 5.0)
    s2 = extrema_over_threshold(y2, 0.3)
    @test length(s2.x) == 2
    @test s2.x_tot[1] == s2.x_tot[2]
    # a flat step on a rising edge is no site, a flat top is one site
    @test length(extrema_over_threshold([0.0, 1.0, 2.0, 2.0, 3.0, 1.0, 0.0], 0.5).x) == 1
    @test length(extrema_over_threshold([0.0, 1.0, 3.0, 3.0, 3.0, 1.0, 0.0], 0.5).x) == 1
    @test_throws ArgumentError extrema_over_threshold(y, -0.3)
    @test_throws ArgumentError extrema_over_threshold(y, NaN)
end

@testset "SiteSolver" begin
    n = 600
    response = bump(n, 300, 12.0)
    solver = SiteSolver(response = response, i_peak = 300, merge_distance = 10.0)
    # two overlapping photons of 1 and 2 units, 15 samples apart, plus a negative site far away
    A_true = [1.0, 2.0, -0.8]
    pos = [200.0, 215.0, 450.0]
    y = sum(A_true[k] .* bump(n, pos[k], 12.0) for k in 1:3)
    sol = solver(y, pos, [5.0, 5.0, 5.0], [1, 1, -1])
    @test sol.A ≈ A_true atol = 1e-6
    @test sol.readoff ≉ A_true                      # the read values carry the neighbour's contribution
    @test sol.pos == pos
    # two same-sign sites closer than merge_distance become one site
    sol2 = solver(y, [200.0, 206.0, 450.0], [5.0, 7.0, 5.0], [1, 1, -1])
    @test length(sol2.A) == 2
    @test sol2.tot[1] == 7.0
    @test 200.0 <= sol2.pos[1] <= 206.0
    # a negative site whose unconstrained solution is positive is set to zero, the rest solved again
    y3 = bump(n, 200, 12.0) .+ 0.5 .* bump(n, 215, 12.0)
    sol3 = solver(y3, [200.0, 215.0], [5.0, 5.0], [1, -1])
    @test sol3.A[2] == 0.0
    @test sol3.A[1] ≈ y3[200] atol = 1e-9
    # an asymmetric response fixes the orientation of the matrix: K[i, k] = response(pos_i - pos_k)
    n4 = 400
    resp = [d < 0 ? exp(-0.5 * (d / 8)^2) : exp(-0.5 * (d / 20)^2) for d in (1:n4) .- 200]
    shifted(k) = [1 <= 200 + (i - k) <= n4 ? resp[200 + (i - k)] : 0.0 for i in 1:n4]
    y4 = 1.0 .* shifted(100) .+ 2.0 .* shifted(130) .- 0.5 .* shifted(170)
    sol4 = SiteSolver(response = resp, i_peak = 200, merge_distance = 5.0)(y4, [100.0, 130.0, 170.0], [5.0, 5.0, 5.0], [1, 1, -1])
    @test sol4.A ≈ [1.0, 2.0, -0.5] atol = 1e-9
    # positions are sorted on output and must lie within the trace
    @test issorted(solver(y, [450.0, 200.0, 215.0], [5.0, 5.0, 5.0], [-1, 1, 1]).pos)
    @test_throws ArgumentError solver(y, [200.0, n + 5.0], [5.0, 5.0], [1, 1])
end
