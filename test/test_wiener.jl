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
