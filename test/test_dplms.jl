using Test
using LegendDSP
using RadiationDetectorDSP
using RadiationDetectorSignals
using Unitful
using IntervalSets
using OffsetArrays
using LinearAlgebra
using Statistics
using Random

@testset "get_wvf_maximum" begin
    Δt = 128.0u"ns"
    t = (0:99) .* Δt
    # parabola with vertex between samples 50 and 51 at t = 50.3 Δt, value 2.0
    y = [2.0 - 0.01 * (i - 50.3)^2 for i in 0:99]
    wvf = RDWaveform(t, y)
    r = get_wvf_maximum(wvf, 40 * Δt, 60 * Δt)
    @test r.max ≈ 2.0 atol = 1e-12
    @test r.t ≈ 50.3 * Δt atol = 1e-9u"ns"

    # plain samples: time is a fractional index (1-based)
    rs = get_wvf_maximum(y, 41, 61)
    @test rs.max ≈ 2.0 atol = 1e-12
    @test rs.t ≈ 51.3 atol = 1e-9

    # offset axes: identical value, index follows the array's own axes
    yo = OffsetArray(y, -10)   # indices -9:90
    ro = get_wvf_maximum(yo, 31, 51)
    @test ro.max ≈ 2.0 atol = 1e-12
    @test ro.t ≈ 41.3 atol = 1e-9

    # a view behaves like the parent
    rv = get_wvf_maximum(view(y, :), 41, 61)
    @test rv.max == rs.max && rv.t == rs.t

    # maximum at the window edge: discrete sample, no interpolation
    re = get_wvf_maximum(wvf, 0 * Δt, 45 * Δt)
    @test re.t == 45 * Δt
    @test re.max == y[46]

    # flat window: no division by zero
    rf = get_wvf_maximum(RDWaveform(t, fill(1.0, 100)), 10 * Δt, 20 * Δt)
    @test rf.max == 1.0 && isfinite(ustrip(rf.t))

    # interval leaves the record on either side
    @test_throws ArgumentError get_wvf_maximum(wvf, -20 * Δt, 10 * Δt)
    @test_throws ArgumentError get_wvf_maximum(wvf, 90 * Δt, 120 * Δt)

    # mismatched time axis and signal axes
    @test_throws DimensionMismatch LegendDSP._get_wvf_maximum_impl(1:10, ones(9), 2:5)
end

@testset "align_waveform" begin
    Δt = 128.0u"ns"
    t = (0:199) .* Δt
    pulse(x) = 1 / (1 + exp(-(x - 100) / 3))        # smooth rise centered at sample 100
    y = [pulse(i) for i in 0:199]
    y_shifted = [pulse(i - 0.37) for i in 0:199]     # same pulse, 0.37 samples later
    wvf = RDWaveform(t, y_shifted)
    a = align_waveform(wvf, 100.37 * Δt, 100 * Δt)
    @test a.time === wvf.time
    @test length(a.signal) == length(y)
    @test maximum(abs.(a.signal[20:180] .- y[20:180])) < 2e-4

    # integer shift of 3 samples to the left, edges repeat the first/last sample
    b = align_waveform(RDWaveform(t, y), 103 * Δt, 100 * Δt)
    @test b.signal[1:190] ≈ y[4:193]
    @test all(b.signal[198:200] .≈ y[200])
    c = align_waveform(RDWaveform(t, y), 97 * Δt, 100 * Δt)
    @test all(c.signal[1:3] .≈ y[1])

    # offset axes are preserved
    yo = OffsetArray(y, -50)
    so = LegendDSP._shift_samples(yo, 0.5)
    @test axes(so) == axes(yo)
    @test so[0] ≈ (LegendDSP._shift_samples(y, 0.5))[50]   # yo[0] is y[50]

    # batch version
    wvfs = ArrayOfRDWaveforms([RDWaveform(t, y_shifted), RDWaveform(t, y)])
    al = align_waveforms(wvfs, [100.37, 100.0] .* Δt, 100 * Δt)
    @test al isa ArrayOfRDWaveforms
    @test maximum(abs.(al.signal[1][20:180] .- al.signal[2][20:180])) < 2e-4
    @test_throws DimensionMismatch align_waveforms(wvfs, [100.0] .* Δt, 100 * Δt)

    # window extraction
    W = aligned_windows(al, 100 * Δt, 10, 20)
    @test size(W) == (31, 2)
    @test W[:, 2] == y[91:121]
    @test_throws ArgumentError aligned_windows(al, 5 * Δt, 10, 20)

    # waveforms without one shared time axis
    wvfs_mixed = ArrayOfRDWaveforms([RDWaveform(t, y), RDWaveform(t .+ Δt, y)])
    @test_throws ArgumentError aligned_windows(wvfs_mixed, 100 * Δt, 10, 20)
end

@testset "dplms_noise_matrix" begin
    rng = MersenneTwister(1)
    Δt = 128.0u"ns"
    t = (0:1023) .* Δt
    n = 20
    # white noise, σ = 3: matrix → 9 I
    wvfs = ArrayOfRDWaveforms([RDWaveform(t, 3 .* randn(rng, 1024)) for _ in 1:400])
    N = dplms_noise_matrix(wvfs, 0.0u"µs" .. 39.0u"µs", n)
    @test size(N) == (n, n)
    @test N ≈ N'
    @test all(abs.(diag(N) .- 9) .< 0.5)
    @test maximum(abs.(N .- Diagonal(diag(N)))) < 0.5

    # AR(1) noise x_i = ρ x_{i-1} + ε: covariance ρ^|i-j| σ²/(1-ρ²)
    ρ = 0.7
    function ar1(rng, len)
        x = zeros(len); x[1] = randn(rng) / sqrt(1 - ρ^2)
        for i in 2:len; x[i] = ρ * x[i-1] + randn(rng); end
        x
    end
    wvfs_ar = ArrayOfRDWaveforms([RDWaveform(t, ar1(rng, 1024)) for _ in 1:400])
    Nar = dplms_noise_matrix(wvfs_ar, 0.0u"µs" .. 39.0u"µs", n)
    expected = [ρ^abs(i - j) / (1 - ρ^2) for i in 1:n, j in 1:n]
    @test maximum(abs.(Nar .- expected)) < 0.15

    # baseline window shorter than the filter
    @test_throws ArgumentError dplms_noise_matrix(wvfs, 0.0u"µs" .. 1.0u"µs", n)
end

@testset "dplms reference and flat-top matrices" begin
    n = 31
    step_pulse(a) = [i < 16 ? 0.0 : a for i in 1:n]           # rise between rows 15 and 16
    windows = hcat(step_pulse(100.0), step_pulse(120.0), step_pulse(80.0))
    r = dplms_reference(windows, 25:31)
    @test r.plateau ≈ 100.0
    @test r.ref[1:15] == zeros(15)
    @test all(r.ref[16:31] .≈ 1.0)
    @test_throws ArgumentError dplms_reference(-windows, 25:31)

    # flat top: a step has slope only at the rise, so F is zero away from it
    ff = 1
    wide = vcat(zeros(1, 3), windows, windows[end:end, :])     # widened by ff on each side
    F = dplms_flattop_matrix(wide, ff)
    @test size(F) == (n, n)
    @test F ≈ F'
    # rows/cols far from the rise are zero; the rise block is positive
    @test all(F[25:31, 25:31] .== 0)
    @test findall(!iszero, F) == [CartesianIndex(15, 15), CartesianIndex(16, 16)]
    @test_throws ArgumentError dplms_flattop_matrix(wide, 0)

    # binned reference with a common plateau
    bins = [1, 2, 2]
    b = dplms_binned_reference(windows, bins, r.plateau)
    @test b.vec[31] ≈ (1.0 + 1.0) / 2                          # bin means 1.0 and 1.0
    @test size(b.mat) == (n, n)
    @test b.mat ≈ (r.ref * r.ref' + r.ref * r.ref') / 2 atol = 1e-12
    @test size(b.refs) == (n, 2)
    @test b.refs[:, 1] ≈ r.ref
    @test_throws ArgumentError dplms_binned_reference(windows, [1, 3, 3], r.plateau)
    @test_throws DimensionMismatch dplms_binned_reference(windows, [1, 2], r.plateau)
end

@testset "dplms_synthesis" begin
    n = 31
    ref = [i < 16 ? 0.0 : 1.0 for i in 1:n]
    σ = 0.01
    noise = σ^2 * Matrix{Float64}(I, n, n)
    flattop = zeros(n, n)
    binned = (mat = zeros(n, n), vec = zeros(n))

    # white noise, no other penalty: matched filter, response 1, noise σ·‖x‖
    s = dplms_synthesis(noise, ref, flattop, binned; nm = 1.0, za = 0.0, ft = 0.0)
    x = reverse(s.coeffs)
    @test dot(x, ref) ≈ 1.0
    @test s.noise_out ≈ σ * norm(x)
    @test all(x[1:15] .≈ 0) && all(x[16:31] .≈ 1 / 16)

    # zero-area penalty drives the coefficient sum to zero
    z = dplms_synthesis(noise, ref, flattop, binned; nm = 1.0, za = 1e6, ft = 0.0)
    @test abs(z.zero_area) < 1e-4
    @test dot(reverse(z.coeffs), ref) ≈ 1.0

    # the ConvolutionFilter reads the window energy at the window end
    Δt = 128.0u"ns"
    t = (0:199) .* Δt
    y = [i < 100 ? 0.0 : 250.0 for i in 0:199]     # step of 250 at sample index 101
    flt = ConvolutionFilter(DirectConvolution(), z.coeffs)
    out = flt(RDWaveform(t, y))
    # window of 31 samples ending at index 101 + 15 = 116 spans the rise like `ref`
    i_end = 116
    @test out.signal[findfirst(==(t[i_end]), out.time)] ≈ 250.0

    # binned term pulls two bins with different plateaus toward the same response
    ref_lo = 0.98 .* ref; ref_lo[10:15] .= 0.02      # lower plateau, slower start
    ref_hi = 1.02 .* ref
    b = (mat = (ref_lo * ref_lo' + ref_hi * ref_hi') / 2, vec = (ref_lo + ref_hi) / 2)
    s0 = dplms_synthesis(noise, ref, flattop, b; nm = 1.0, za = 0.0, ft = 0.0, ct = 0.0)
    s1 = dplms_synthesis(noise, ref, flattop, b; nm = 1.0, za = 0.0, ft = 0.0, ct = 1e3)
    gap(s) = abs(dot(reverse(s.coeffs), ref_hi) - dot(reverse(s.coeffs), ref_lo))
    @test gap(s1) < gap(s0)

    # size mismatch and singular system fail fast
    @test_throws DimensionMismatch dplms_synthesis(noise[1:10, 1:10], ref, flattop, binned; nm = 1.0, za = 0.0, ft = 0.0)
    @test_throws ArgumentError dplms_synthesis(zeros(n, n), ref, flattop, binned; nm = 1.0, za = 0.0, ft = 0.0)
end

@isdefined(make_fake_config) || include("fake_config.jl")

# Presummed-like pulse: 1024 samples at 128 ns, logistic rise of width `w` samples centered
# at 50 µs (so half of the current integral is reached there), exponential tail with decay `τ`,
# white baseline noise σ, DC offset 1000.
function make_dplms_waveform(rng; w = 2.0, amplitude = 1e4, τ = 400.0u"µs", σ = 5.0)
    Δt = 128.0u"ns"
    n = 1024
    t = (0:n-1) .* Δt
    i_c = 1 + round(Int, 50.0u"µs" / Δt)
    τs = ustrip(NoUnits, τ / Δt)
    y = σ .* randn(rng, n)
    for i in eachindex(y)
        y[i] += amplitude / (1 + exp(-(i - i_c) / w)) * exp(-max(i - i_c, 0) / τs)
    end
    RDWaveform(t, y .+ 1000.0)
end

@testset "dsp_dplms_optimization" begin
    rng = MersenneTwister(7)
    widths = [isodd(k) ? 2.0 : 3.0 for k in 1:300]
    wvfs = ArrayOfRDWaveforms([make_dplms_waveform(rng; w = w) for w in widths])
    config = make_fake_config()
    τ = 400.0u"µs"
    res = dsp_dplms_optimization(wvfs, config, τ, widths; n_pre = 40, n_post = 40, n_bins = 2)
    g = length(res.grid)
    @test size(res.e_grid) == (g, 300)
    @test size(res.emax_grid) == (g, 300)
    @test length(res.coeffs) == g && all(length.(res.coeffs) .== 81)
    @test res.t_pickoff == 40 * 128.0u"ns"
    @test all(isfinite, res.e_grid)
    # every grid point reconstructs the amplitude to better than 1 percent on average
    @test all(abs.(mean(res.e_grid, dims = 2) .- 1e4) .< 100)
    # predicted noise is small compared with the amplitude and positive
    @test all(0 .< res.noise_out .< 0.01)
    @test length(res.ref) == 81
    @test res.ref[end] ≈ 1.0 atol = 0.05
    @test size(res.ref_bins, 2) == 2
    @test length(res.t50) == 300
    # t50 lands within one sample of the rise center for every event
    @test all(abs.(res.t50 .- 50.0u"µs") .<= 128.0u"ns")
    # t50 is the t50_pre of dsp_icpc_compressed
    bl = signalstats.(wvfs, leftendpoint(config.bl_window), rightendpoint(config.bl_window))
    wvfs_pz = InvCRFilter(τ).(shift_waveform.(wvfs, -bl.mean))
    e_101010 = extremestats.(TrapezoidalChargeFilter(10u"µs", 10u"µs").(wvfs_pz)).max
    @test res.t50 == get_threshold(wvfs_pz, e_101010 .* 0.5; mintot = config.kwargs_pars.tx_mintot)
    # widenings beyond n_post leave no flat region in the 81-sample window
    @test length(res.flat_top) == length(res.shaping_time) == g
    @test res.flat_top[1] > 0u"ns" && all(>=(0u"ns"), res.flat_top) && all(>(0u"ns"), res.shaping_time)
    @test length(res.superpulse_bins) == 2
    # superpulse and per-grid-point response
    @test length(res.response) == g
    @test res.response[1] isa RDWaveform
    @test size(res.noise) == (81, 81) && all(size(F) == (81, 81) for F in values(res.flattops))
    @test res.noise ≈ res.noise'
    @test all(F ≈ F' for F in values(res.flattops))
    # the w = 3 logistic reaches 0.999 of its plateau 21 samples after its center
    @test res.n_rise in 20:22
    @test [w.ff for w in res.grid] == res.n_rise .* (1:4)
    r_noza = dsp_dplms_optimization(wvfs, config, τ, widths; n_pre = 40, n_post = 40, n_bins = 2,
        grid = n -> dplms_default_grid(n; za = 0.0, ff_factors = (1, 1.5)))
    @test all(w.za == 0 for w in r_noza.grid) && [w.ff for w in r_noza.grid] == round.(Int, r_noza.n_rise .* [1, 1.5])
    @test_throws "never reaches" dsp_dplms_optimization(wvfs, config, τ, widths; n_pre = 40, n_post = 40, n_bins = 2, ff_level = 2.0)
    se = SignalEstimator(PolynomialDNI(config.kwargs_pars.sig_interpolation_order, config.kwargs_pars.sig_interpolation_length))
    t = first(wvfs.time)
    t_target = t[LegendDSP._sample_index(t, median(res.t50))]
    @test se(res.response[1], t_target + res.t_pickoff) ≈ 1.0 atol = 0.05
    # noise_out predicts the output spread of the filter on pure noise of the training σ
    r0 = dsp_dplms_optimization(wvfs, config, τ, widths; n_pre = 40, n_post = 40, n_bins = 2,
        grid = [(nm = 1.0, za = 0.0, ft = 0.0, ct = 0.0, ff = 1)])
    noise_wvfs = ArrayOfRDWaveforms([RDWaveform((0:1023) .* 128.0u"ns", 5.0 .* randn(rng, 1024)) for _ in 1:200])
    flt_noise = ConvolutionFilter(DirectConvolution(), r0.coeffs[1]).(noise_wvfs)
    # dplms_noise_matrix underestimates the lowest-frequency noise power because the training
    # waveforms' baseline mean was subtracted over the same window; the za = ft = 0 filter has
    # a large DC component (zero_area ≈ 0.6-0.8), so the prediction comes out about 11% low
    @test std(reduce(vcat, flt_noise.signal)) ≈ r0.noise_out[1] * r0.plateau rtol = 0.15
    # a ct > 0 grid point runs and stays finite
    res_ct = dsp_dplms_optimization(wvfs, config, τ, widths; n_pre = 40, n_post = 40, n_bins = 2,
        grid = [(nm = 1.0, za = 0.0, ft = 0.0, ct = 1e2, ff = 1)])
    @test all(isfinite, res_ct.e_grid)
    # drift vector of the wrong length fails
    @test_throws DimensionMismatch dsp_dplms_optimization(wvfs, config, τ, [1.0, 2.0]; n_pre = 40, n_post = 40)
    @test_throws "must be finite" dsp_dplms_optimization(wvfs, config, τ, [NaN; widths[2:end]]; n_pre = 40, n_post = 40)
    # a waveform without a rising edge has no t50 inside the current window
    flat_wvf = RDWaveform(first(wvfs.time), fill(1000.0, 1024))
    @test_throws "for 1 of 300 waveforms" dsp_dplms_optimization(ArrayOfRDWaveforms([flat_wvf; wvfs[2:end]]), config, τ, widths; n_pre = 40, n_post = 40)
end

@testset "dsp_dplms_optimization default window" begin
    rng = MersenneTwister(7)
    widths = [isodd(k) ? 2.0 : 3.0 for k in 1:100]
    wvfs = ArrayOfRDWaveforms([make_dplms_waveform(rng; w = w) for w in widths])
    config = make_fake_config()
    Δt = 128.0u"ns"
    n_bl = round(Int, config.bl_window.right / Δt) - round(Int, config.bl_window.left / Δt) + 1
    grid = [(nm = 1.0, za = 1e3, ft = 1.0, ct = 0.0, ff = ff) for ff in (1, 3)]
    res = dsp_dplms_optimization(wvfs, config, 400.0u"µs", widths; n_bins = 2, grid)
    @test res.n_pre == res.n_post == (n_bl - 8) ÷ 2
    @test all(length.(res.coeffs) .== 2 * res.n_pre + 1)
    @test all(isfinite, res.e_grid)
    @test res.t_pickoff == res.n_post * Δt
    # wider flat-top terms give longer flat tops and shorter shaping times
    r_ff = dsp_dplms_optimization(wvfs, config, 400.0u"µs", widths; n_bins = 2)
    @test all(>(0u"ns"), r_ff.flat_top) && issorted(r_ff.flat_top) && issorted(r_ff.shaping_time; rev = true)
    @test_throws ArgumentError dsp_dplms_optimization(wvfs, config, 400.0u"µs", widths; n_pre = 200, n_post = 200, grid)
end

@testset "dsp_dplms_energy" begin
    rng = MersenneTwister(3)
    wvfs = ArrayOfRDWaveforms([make_dplms_waveform(rng) for _ in 1:50])
    config = make_fake_config()
    # a flat-topped filter whose read-out time t50 + t_pickoff lies after the end of the current window
    res = dsp_dplms_optimization(wvfs, config, 400.0u"µs", collect(1.0:50.0); n_pre = 100, n_post = 140, n_bins = 1,
        grid = n -> dplms_default_grid(n)[1:1])
    @test all(>(rightendpoint(config.current_window)), res.t50 .+ res.t_pickoff)
    bl = signalstats.(wvfs, 0.0u"µs", 39.0u"µs")
    wvfs_pz = InvCRFilter(400.0u"µs").(shift_waveform.(wvfs, -bl.mean))
    se = SignalEstimator(PolynomialDNI(config.kwargs_pars.sig_interpolation_order, config.kwargs_pars.sig_interpolation_length))
    e = dsp_dplms_energy(wvfs_pz, res.coeffs[1], res.t50, res.t_pickoff, config.current_window, se)
    @test e.e ≈ vec(res.e_grid[1, :])
    @test e.emax ≈ vec(res.emax_grid[1, :])
    # the maximum is the flat top around the read-out time
    @test all(e.emax .>= 0.999 .* e.e) && median(e.emax ./ e.e) < 1.001
end

@testset "get_dplms_pars" begin
    config = make_fake_config()
    Δt = 128.0u"ns"
    stored = get_dplms_pars(PropDict(:dplms => PropDict(:coefficients => [0.1, 0.2, 0.7], :t_pickoff => 1.0u"µs")), config, Δt)
    @test stored.coefficients == [0.1, 0.2, 0.7] && stored.t_pickoff == 1.0u"µs"
    # the default filter is the trapezoid of flt_defaults.dplms read out like e_trap
    default = get_dplms_pars(PropDict(), config, Δt)
    @test default.t_pickoff == 6.5u"µs"
    rng = MersenneTwister(5)
    wvfs = ArrayOfRDWaveforms([make_dplms_waveform(rng) for _ in 1:20])
    wvfs_pz = InvCRFilter(400.0u"µs").(shift_waveform.(wvfs, -1000.0))
    t50 = fill(50.0u"µs", 20)
    se = SignalEstimator(PolynomialDNI(config.kwargs_pars.sig_interpolation_order, config.kwargs_pars.sig_interpolation_length))
    @test dsp_dplms_energy(wvfs_pz, default.coefficients, t50, default.t_pickoff, config.current_window, se).e ≈
        se.(TrapezoidalChargeFilter(5.0u"µs", 3.0u"µs").(wvfs_pz), t50 .+ 6.5u"µs") rtol = 1e-6
    pd_nodefault = make_fake_config_pd()
    delete!(pd_nodefault.flt_defaults, :dplms)
    @test_throws "flt_defaults.dplms" get_dplms_pars(PropDict(), DSPConfig(pd_nodefault), Δt)
    @test_throws "no coefficients" get_dplms_pars(PropDict(:dplms => PropDict(:coefficients => Float64[], :t_pickoff => 1.0u"µs")), config, Δt)
    @test_throws "no t_pickoff" get_dplms_pars(PropDict(:dplms => PropDict(:coefficients => [1.0])), config, Δt)
end
