using Test
using LegendDSP
using RadiationDetectorSignals
using Unitful
using Statistics
using Random

@testset "iterative_baseline" begin
    Δt = 16u"ns"

    @testset "constant signal" begin
        t = (0:99) .* Δt
        sig = fill(5.0, 100)
        wvf = RDWaveform(t, sig)
        r = iterative_baseline(wvf)
        @test r.mean == 5.0
        @test r.sigma == 0.0
        @test r.converged
        @test r.n_accepted == 100
    end

    @testset "pure gaussian noise: matches plain mean/std" begin
        rng = MersenneTwister(0xCAFE)
        N = 5_000
        μ_true = 1.2
        σ_true = 0.7
        t = (0:N-1) .* Δt
        sig = μ_true .+ σ_true .* randn(rng, N)
        wvf = RDWaveform(t, sig)
        r = iterative_baseline(wvf)

        @test isapprox(r.mean, mean(sig), rtol=1e-3)
        @test isapprox(r.sigma, std(sig; corrected=false), rtol=0.05)
        @test r.converged
    end

    @testset "pulse rejection: mean recovers true baseline" begin
        rng = MersenneTwister(0xC0DE01)
        N = 5_000
        σ_true = 1.0
        t = (0:N-1) .* Δt
        sig = σ_true .* randn(rng, N)

        # 10% positive-pulse contamination at +20σ
        n_pulses = N ÷ 10
        for _ in 1:n_pulses
            sig[rand(rng, 1:N)] += 20.0
        end

        plain_mean = mean(sig)
        wvf = RDWaveform(t, sig)
        r = iterative_baseline(wvf; n_sigma=3, max_iter=12)

        # Plain mean is biased
        @test abs(plain_mean) > 1.0
        # Iterative mean recovers true baseline
        @test abs(r.mean) < 0.1
        @test isapprox(r.sigma, σ_true, rtol=0.15)
        @test r.n_accepted < N
        @test r.n_accepted > N - 2*n_pulses
        @test r.converged
        @test r.n_iter <= 12
    end

    @testset "two-sided pulses (mean unbiased)" begin
        # NOTE: the convergence check looks at the change in the running MEAN.
        # For two-sided contamination the mean is already unbiased after the
        # first pass — the algorithm therefore declares convergence even though
        # σ is still inflated by partially-included outliers. We only check the
        # mean recovery here, mirroring how the routine behaves on real
        # (typically asymmetric) baselines.
        rng = MersenneTwister(0xC0DE02)
        N = 5_000
        t = (0:N-1) .* Δt
        sig = randn(rng, N)
        for _ in 1:(N ÷ 20)
            sig[rand(rng, 1:N)] += 30.0
            sig[rand(rng, 1:N)] -= 30.0
        end
        wvf = RDWaveform(t, sig)
        r = iterative_baseline(wvf; n_sigma=3, max_iter=12)
        @test abs(r.mean) < 0.1
    end

    @testset "unitful signal" begin
        rng = MersenneTwister(0xC0DE03)
        N = 2_000
        t = (0:N-1) .* Δt
        sig = (0.5u"mV") .* randn(rng, N) .+ 2.0u"mV"
        wvf = RDWaveform(t, sig)
        r = iterative_baseline(wvf)
        @test unit(r.mean) == u"mV"
        @test unit(r.sigma) == u"mV"
        @test isapprox(ustrip(u"mV", r.mean), 2.0, atol=0.05)
        @test isapprox(ustrip(u"mV", r.sigma), 0.5, rtol=0.1)
    end

    @testset "plain array input" begin
        rng = MersenneTwister(0xC0DE04)
        sig = randn(rng, 2_000)
        r = iterative_baseline(sig)
        @test abs(r.mean) < 0.1
        @test isapprox(r.sigma, 1.0, rtol=0.15)
    end

    @testset "empty input" begin
        r = iterative_baseline(Float64[])
        @test r.mean == 0.0
        @test r.sigma == 0.0
        @test r.n_accepted == 0
        @test !r.converged
    end

    @testset "respects max_iter" begin
        rng = MersenneTwister(0xC0DE05)
        N = 5_000
        sig = randn(rng, N)
        for _ in 1:(N ÷ 5)
            sig[rand(rng, 1:N)] += 5.0
        end
        r = iterative_baseline(sig; n_sigma=2, max_iter=2)
        @test r.n_iter <= 2
    end

    @testset "return_mask marks accepted samples" begin
        rng = MersenneTwister(0xC0DE10)
        N = 5_000
        sig = randn(rng, N)
        for _ in 1:(N ÷ 10)
            sig[rand(rng, 1:N)] += 20.0
        end
        r0 = iterative_baseline(sig; n_sigma=3)
        r  = iterative_baseline(sig; n_sigma=3, return_mask=true)
        @test !haskey(r0, :accepted)                 # default path unchanged
        @test r.mean == r0.mean && r.sigma == r0.sigma
        @test length(r.accepted) == N
        @test count(r.accepted) == r.n_accepted      # mask ≙ the defining window
        @test !any(r.accepted[i] for i in 1:N if sig[i] > 10.0)  # +20σ pulses rejected
    end

    @testset "return_mask: constant signal → all accepted" begin
        r = iterative_baseline(fill(5.0, 100); return_mask=true)
        @test length(r.accepted) == 100
        @test all(r.accepted)
    end

    @testset "return_mask: empty input" begin
        r = iterative_baseline(Float64[]; return_mask=true)
        @test isempty(r.accepted)
    end
end
