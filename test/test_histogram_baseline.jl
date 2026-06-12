using Test
using LegendDSP
using RadiationDetectorSignals
using Unitful
using Statistics
using Random

@testset "histogram_baseline" begin
    Δt = 16u"ns"

    @testset "constant signal" begin
        t = (0:99) .* Δt
        sig = fill(5.0, 100)
        wvf = RDWaveform(t, sig)
        r = histogram_baseline(wvf)
        @test r.mode ≈ 5.0 atol=1e-10
        @test r.sigma == 0.0
        @test r.hwhm_left == 0.0
        @test r.hwhm_right == 0.0
        @test r.n_used == 100
        @test r.mode_count == 100
    end

    @testset "gaussian noise: mode and HWHM recovery" begin
        # Fresh RNG so the test is independent of test ordering.
        rng = MersenneTwister(0xBA5E11)
        N = 50_000
        μ_true = 2.5
        σ_true = 1.0
        t = (0:N-1) .* Δt
        sig = μ_true .+ σ_true .* randn(rng, N)
        wvf = RDWaveform(t, sig)
        r = histogram_baseline(wvf; nbins=200)

        # Smoothed mode + parabolic interpolation typically gets within
        # ~σ/15 of the true mean for 50k samples.
        @test isapprox(r.mode, μ_true, atol=0.1)

        k = sqrt(2 * log(2))
        @test isapprox(r.hwhm_left,  k * σ_true, rtol=0.10)
        @test isapprox(r.hwhm_right, k * σ_true, rtol=0.10)
        @test isapprox(r.sigma, σ_true, rtol=0.10)

        # Symmetric distribution → no significant asymmetry
        @test abs(r.hwhm_right - r.hwhm_left) / r.sigma < 0.2
    end

    @testset "asymmetry detection (pulse contamination)" begin
        # Heavy near-baseline contamination on the positive side. Pulses
        # land mostly within ±2σ of the noise so they actually distort the
        # half-max region (the part HWHM probes).
        rng = MersenneTwister(0xCAFE01)
        N = 50_000
        σ_true = 1.0
        t = (0:N-1) .* Δt
        sig = σ_true .* randn(rng, N)
        # 30 % of samples acquire a positive offset of ~1–2σ — strong enough
        # to widen the right flank of the histogram measurably.
        n_pulses = 15_000
        for _ in 1:n_pulses
            i = rand(rng, 1:N)
            sig[i] += 1.5 + 0.5 * randn(rng)
        end
        wvf = RDWaveform(t, sig)
        r = histogram_baseline(wvf; nbins=200)

        # Mode still near zero — the unaffected 70 % bulk still dominates.
        @test abs(r.mode) < 0.2

        # Right HWHM clearly larger than left HWHM.
        @test r.hwhm_right > r.hwhm_left * 1.10
    end

    @testset "range bounds drop pulses" begin
        rng = MersenneTwister(0xCAFE02)
        N = 20_000
        t = (0:N-1) .* Δt
        sig = randn(rng, N)
        # large outliers far outside the analysis window
        sig[1:50] .= 100.0
        wvf = RDWaveform(t, sig)
        r = histogram_baseline(wvf; nbins=100, min=-5.0, max=5.0)
        @test r.n_used == N - 50
        @test abs(r.mode) < 0.2          # 1-bin-width tolerance
        @test isapprox(r.sigma, 1.0, rtol=0.15)
    end

    @testset "unitful signal" begin
        rng = MersenneTwister(0xCAFE03)
        N = 10_000
        t = (0:N-1) .* Δt
        sig = (1.0u"mV") .* randn(rng, N) .+ 3.0u"mV"
        wvf = RDWaveform(t, sig)
        r = histogram_baseline(wvf; nbins=100)
        @test unit(r.mode) == u"mV"
        @test unit(r.sigma) == u"mV"
        @test isapprox(ustrip(u"mV", r.mode), 3.0, atol=0.15)
        @test isapprox(ustrip(u"mV", r.sigma), 1.0, rtol=0.15)
    end

    @testset "plain array input" begin
        rng = MersenneTwister(0xCAFE04)
        sig = randn(rng, 20_000)
        r = histogram_baseline(sig; nbins=100)
        @test abs(r.mode) < 0.2
        @test isapprox(r.sigma, 1.0, rtol=0.15)
    end

    @testset "smoothing disabled" begin
        # When smooth=1 the algorithm uses raw histogram bins: mode is
        # noisier but still recovers true mean within a few bin widths.
        rng = MersenneTwister(0xCAFE05)
        N = 50_000
        sig = randn(rng, N)
        r = histogram_baseline(sig; nbins=200, smooth=1)
        @test abs(r.mode) < 0.3
        @test isapprox(r.sigma, 1.0, rtol=0.20)
    end

    @testset "degenerate: empty range" begin
        sig = fill(5.0, 100)
        r = histogram_baseline(sig; nbins=100, min=10.0, max=20.0)
        @test r.n_used == 0
        @test r.mode == 0.0
        @test r.sigma == 0.0
    end

    @testset "robust window rejects extreme outliers" begin
        # Default behaviour: range = median ± n_sigma·MAD, so saturated samples
        # far from the baseline neither move the centre nor inflate the scale —
        # the bins stay narrow and the baseline is recovered without any min/max.
        rng = MersenneTwister(0xCAFE06)
        N = 20_000
        sig = randn(rng, N)
        sig[1:10] .= 16_000.0
        r = histogram_baseline(sig; nbins=100)
        @test abs(r.mode) < 0.2
        @test isapprox(r.sigma, 1.0, rtol=0.15)
    end

    @testset "fixed bin width" begin
        rng = MersenneTwister(0xCAFE07)
        N = 50_000
        sig = 2.0 .+ randn(rng, N)
        sig[1] = 1_000.0   # huge range, but the bin width stays fixed
        r = histogram_baseline(sig; bin_width=0.05)
        @test isapprox(r.mode, 2.0, atol=0.1)
        @test isapprox(r.sigma, 1.0, rtol=0.15)
    end

    @testset "interpolation disabled" begin
        rng = MersenneTwister(0xCAFE08)
        N = 50_000
        sig = randn(rng, N)
        r = histogram_baseline(sig; nbins=200, interpolation=false)
        @test abs(r.mode) < 0.3   # nearest-bin mode: coarser by design (no parabola)
        @test isapprox(r.sigma, 1.0, rtol=0.20)
        @test r.hwhm_left > 0
        @test r.hwhm_right > 0
    end

    @testset "return_hist exposes the binned histogram" begin
        rng = MersenneTwister(0xBA5E22)
        N = 50_000
        sig = 2.5 .+ randn(rng, N)
        r0 = histogram_baseline(sig; nbins=200)
        r  = histogram_baseline(sig; nbins=200, return_hist=true)
        @test !haskey(r0, :hist)                       # default path unchanged
        @test r.mode == r0.mode && r.sigma == r0.sigma
        @test length(r.hist.counts) == 200
        @test length(r.hist.smoothed) == 200
        @test length(r.hist.centers) == 200
        @test 0 < sum(r.hist.counts) <= r.n_used       # robust window drops tail samples
        @test argmax(r.hist.smoothed) == r.hist.mode_bin
        @test r.hist.half_max ≈ r.hist.smoothed[r.hist.mode_bin] / 2
        @test r.hist.lo < r.mode < r.hist.hi
    end

    @testset "return_hist: fixed bin width" begin
        rng = MersenneTwister(0xBA5E23)
        sig = 2.0 .+ randn(rng, 50_000)
        r = histogram_baseline(sig; bin_width=0.05, return_hist=true)
        @test r.hist.bin_width == 0.05
        @test isapprox(sum(r.hist.counts), r.n_used)   # fixed width: all in-range used
    end

    @testset "return_hist: degenerate empty histogram" begin
        r = histogram_baseline(fill(5.0, 100); min=10.0, max=20.0, return_hist=true)
        @test r.n_used == 0
        @test isempty(r.hist.counts)
        @test isempty(r.hist.centers)
    end
end
