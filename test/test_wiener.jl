using Test
using LegendDSP
using Unitful
using Statistics
using Random

@testset "wiener_template" begin
    rng = MersenneTwister(2)
    n = 200
    shape = [i < 40 ? 0.0 : (1 - exp(-(i - 40) / 2)) * exp(-(i - 40) / 30) for i in 1:n]
    shape ./= maximum(shape)
    heights = 10 .+ 40 .* rand(rng, 500)
    windows = shape .* heights' .+ 0.5 .* randn(rng, n, 500)

    # amplitudes on another scale than the pulse height: seven units of height per unit of amplitude
    r = wiener_template(windows, heights ./ 7)
    @test maximum(r.template) == 1
    @test r.peak ≈ 7 rtol = 0.01
    @test maximum(abs.(r.template .- shape)) < 0.01

    # one window in ten holds a second pulse: the median keeps the shape
    piled = copy(windows)
    piled[120:end, 1:10:end] .+= 30
    @test maximum(abs.(wiener_template(piled, heights ./ 7).template .- shape)) < 0.02

    # integer windows
    v = round.(Int, 100 .* shape)
    @test wiener_template(v .* [1 2 3], [1, 2, 3]).template ≈ v ./ 100

    @test_throws DimensionMismatch wiener_template(windows, heights[1:end-1])
    @test_throws ArgumentError wiener_template(windows, vcat(0.0, heights[2:end]))
    @test_throws ArgumentError wiener_template(windows, vcat(NaN, heights[2:end]))
    @test_throws ArgumentError wiener_template(-(shape .* heights'), heights)
end

@testset "wiener_noise_spectrum" begin
    rng = MersenneTwister(1)
    Δt = 16.0u"ns"
    n = 512

    # white noise, σ = 3: flat spectrum 2 σ² Δt, variance 9
    traces = 3 .* randn(rng, n, 2000)
    s = wiener_noise_spectrum(traces, Δt)
    @test length(s.freq) == length(s.psd) == n ÷ 2 + 1
    @test first(s.freq) == 0.0u"ns^-1"
    @test last(s.freq) ≈ 1 / (2Δt)
    @test unit(step(s.freq)) == u"ns^-1" && unit(first(s.psd)) == u"ns"
    @test sum(s.psd) * step(s.freq) ≈ 9 rtol = 0.02
    @test mean(s.psd[3:end-2]) ≈ 2 * 9 * Δt rtol = 0.02
    @test all(isapprox.(s.psd[3:end-2], 2 * 9 * Δt; rtol = 0.2))

    # moving-average noise x_i = ε_i + a ε_{i-1}: spectrum 2 Δt (1 + a² + 2 a cos(2π f Δt))
    a = 0.8
    e = randn(rng, n + 1, 2000)
    s_ma = wiener_noise_spectrum(e[2:end, :] .+ a .* e[1:end-1, :], Δt)
    expected = 2Δt .* (1 + a^2 .+ 2a .* cospi.(2 .* s_ma.freq .* Δt))
    @test all(isapprox.(s_ma.psd[3:end-2], expected[3:end-2]; rtol = 0.2))

    # one trace in twenty holds a burst: the median hardly moves
    bursts = copy(traces)
    bursts[100:140, 1:20:end] .+= 200
    @test mean(wiener_noise_spectrum(bursts, Δt).psd[3:end-2]) ≈ 2 * 9 * Δt rtol = 0.1

    # a constant and an alternating trace of 8 samples: Hann window with sum 7/2 and squared sum 21/8,
    # the first and the last value are not doubled
    @test wiener_noise_spectrum(ones(8, 1), 1.0).psd[1] ≈ 3.5^2 / (2.625 * log(2))
    @test wiener_noise_spectrum(reshape([(-1.0)^k for k in 0:7], 8, 1), 1.0).psd[end] ≈ 3.5^2 / (2.625 * log(2))

    # plain numbers and the element type of the traces
    s32 = wiener_noise_spectrum(Float32.(traces), 1.0f0)
    @test eltype(s32.psd) == Float32
    @test sum(s32.psd) * step(s32.freq) ≈ 9 rtol = 0.02

    @test_throws ArgumentError wiener_noise_spectrum(traces[1:end-1, :], Δt)
end

@testset "wiener_synthesis" begin
    Δt = 16.0u"ns"
    n, m, n_pre = 512, 200, 50
    σ = 3.0
    # white noise of standard deviation σ, as wiener_noise_spectrum gives it for many traces: N = n σ² at every frequency
    noise = (freq = (0:n÷2) ./ (n * Δt), psd = [k in (0, n ÷ 2) ? σ^2 * Δt : 2 * σ^2 * Δt for k in 0:n÷2])
    lag = (1:m) .- (n_pre + 1)
    target = exp.(-0.5 .* (lag ./ 3) .^ 2)
    target ./= sum(target)

    # a template of one sample: the kernel is the target, scaled by 1 / (1 + λ n σ²)
    spike = [l == 0 ? 1.0 : 0.0 for l in lag]
    λ = 1e-3
    s = wiener_synthesis(noise, spike, target, n_pre; λ, n_half = 40)
    scale = 1 / (1 + λ * n * σ^2)
    @test length(s.coeffs) == 81
    @test s.coeffs ≈ scale .* target[n_pre+1-40:n_pre+1+40]
    @test s.response ≈ scale .* target
    @test s.noise_out ≈ scale * σ * sqrt(sum(abs2, target))
    @test s.spectrum.freq == noise.freq
    @test s.spectrum.template ≈ ones(n ÷ 2 + 1)
    @test s.spectrum.noise ≈ fill(n * σ^2, n ÷ 2 + 1)
    @test first(s.spectrum.target) ≈ 1
    @test s.spectrum.gain ≈ scale .* sqrt.(s.spectrum.target)

    # a pulse with a fast rise and a slow decay: without noise weight the response is the target
    pulse = [l < 0 ? 0.0 : (1 - exp(-(l + 1) / 2)) * exp(-l / 25) for l in lag]
    pulse ./= maximum(pulse)
    full = wiener_synthesis(noise, pulse, target, n_pre; λ = 1e-12, n_half = n ÷ 2 - 1)
    @test maximum(abs.(full.response .- target)) < 1e-6
    @test sum(full.response) ≈ 1

    # the response is the template convolved with the coefficients, lag 0 in the middle
    r = wiener_synthesis(noise, pulse, target, n_pre; λ = 1e-4, n_half = 100)
    conv = [sum(r.coeffs[101+l] * pulse[i-l] for l in -100:100 if 1 <= i - l <= m) for i in 1:m]
    @test r.response ≈ conv atol = 1e-10

    # more noise weight: less noise after the kernel, a response farther from the target
    weak, strong = wiener_synthesis(noise, pulse, target, n_pre; λ = 1e-6, n_half = 100), wiener_synthesis(noise, pulse, target, n_pre; λ = 1e-2, n_half = 100)
    @test strong.noise_out < weak.noise_out
    @test maximum(abs.(strong.response .- target)) > maximum(abs.(weak.response .- target))

    # the predicted noise is that of noise traces through the coefficients, with the spectrum measured from them
    rng = MersenneTwister(3)
    traces = σ .* randn(rng, n, 2000)
    measured = wiener_synthesis(wiener_noise_spectrum(traces, Δt), pulse, target, n_pre; λ = 1e-4, n_half = 100)
    @test measured.coeffs ≈ r.coeffs rtol = 0.05
    filtered = [sum(measured.coeffs[101+l] * traces[i-l, j] for l in -100:100) for i in 101:n-100, j in 1:200]
    @test std(filtered) ≈ measured.noise_out rtol = 0.05

    @test_throws DimensionMismatch wiener_synthesis(noise, pulse, target[1:end-1], n_pre; λ, n_half = 40)
    @test_throws ArgumentError wiener_synthesis(noise, pulse, target, m; λ, n_half = 40)
    @test_throws ArgumentError wiener_synthesis(noise, pulse, target, n_pre; λ, n_half = n ÷ 2)
    @test_throws ArgumentError wiener_synthesis(noise, vcat(pulse, zeros(n)), vcat(target, zeros(n)), n_pre; λ, n_half = 40)
    @test_throws ArgumentError wiener_synthesis(noise, zeros(m), target, n_pre; λ = 0, n_half = 40)
end
