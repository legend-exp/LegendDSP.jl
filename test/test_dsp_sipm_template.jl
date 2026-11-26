using Test
using LegendDSP
using RadiationDetectorSignals
using Unitful
using TypedTables
using PropDicts
using ArraysOfArrays
using Random
using Statistics

@testset "dsp_sipm_pulse_selection_compressed" begin
    n = 6250
    Δt = 16.0u"ns"
    t = range(0.0u"ns", step = Δt, length = n)
    rng = MersenneTwister(4)
    pulse_shape(i0, a) = [i < i0 ? 0.0 : a * (1 - exp(-(i - i0) / 3)) * exp(-(i - i0) / 300) for i in 1:n]
    make_wf(pulses...) = 3750.0 .+ 0.3 .* randn(rng, n) .+ reduce(.+, [pulse_shape(i0, a) for (i0, a) in pulses]; init = zeros(n))
    tail(a) = [a * exp(-i / 1875) for i in 1:n]
    i_mid = round(Int, 48.0u"µs" / Δt)
    signals = vcat(
        [make_wf((i_mid, 20.0)) for _ in 1:10],             #  1-10: one pulse
        [make_wf((i_mid, 40.0))],                           # 11: one pulse of twice the size
        [make_wf((i_mid, 20.0)) .+ tail(1.5)],              # 12: one pulse on the slow tail of an earlier one
        [make_wf((i_mid, 20.0), (i_mid + 200, 20.0))],      # 13: two pulses
        [make_wf((round(Int, 10.0u"µs" / Δt), 20.0))],      # 14: one pulse early in the trace
        [make_wf((round(Int, 80.0u"µs" / Δt), 20.0))],      # 15: one pulse with the window beyond the trace end
        [make_wf() for _ in 1:10],                          # 16-25: noise
        [make_wf() .+ tail(1.0)],                           # 26: noise on the slow tail of an earlier pulse
        [make_wf((i_mid, 2.0))],                            # 27: one pulse above the trigger threshold, below the pulse threshold
    )
    wvfs = ArrayOfRDWaveforms((fill(t, length(signals)), VectorOfSimilarVectors(reduce(hcat, signals))))
    config = PropDict(
        :trigger => PropDict(:mwa_window => 160.0u"ns", :passes => 4, :n_σ_threshold => 4.0, :n_σ_valley => 4.0, :min_tot => 80.0u"ns",
                             :n_σ_pulse => 40.0, :min_tot_pulse => 64.0u"ns",
                             :min_threshold => -2.0, :max_threshold => 2.0),
        :window => [-5.0u"µs", 30.0u"µs"], :baseline_window => [-5.0u"µs", -1.0u"µs"],
        :n_σ_baseline => 3.0, :flat_window => 10.0u"µs", :n_σ_noise => 8.0, :max_noise_wvfs => 5)
    res = dsp_sipm_pulse_selection_compressed(wvfs, config)

    @test res.n_trig == vcat(fill(1, 12), 2, 1, 1, fill(0, 11), 1)
    @test res.is_single == vcat(fill(true, 12), false, true, false, fill(false, 12))
    @test res.is_noise == vcat(fill(false, 15), fill(true, 11), false)
    @test length.(res.pulse) == 2188 .* res.is_single
    @test length.(res.noise) == vcat(fill(0, 15), fill(n, 5), fill(0, 7))

    # time from the start of the trace to the trigger
    @test res.pretrigger[1] ≈ 48.0u"µs" atol = 0.1u"µs"
    @test res.pretrigger[14] ≈ 10.0u"µs" atol = 0.1u"µs"
    @test isnan(res.pretrigger[16])

    # flatness: near zero on a flat trace, large on the slow tail of an earlier pulse, undefined for all other waveforms
    @test all(abs.(res.flatness[[1:11; 16:25]]) .< 0.15)
    @test res.flatness[12] > 0.5 && res.flatness[26] > 0.5
    @test all(isnan, res.flatness[[13, 15, 27]])
    @test 4 < res.trig_amp[27] / res.threshold[27] < 40 < res.trig_amp[1] / res.threshold[1]

    # the trigger sits on the rising edge, and its amplitude scales with the pulse
    @test res.trig_t[1] ≈ 48.0u"µs" atol = 0.1u"µs"
    @test res.trig_amp[11] / res.trig_amp[1] ≈ 2.0 atol = 0.1
    @test isnan(res.trig_t[16]) && isnan(res.trig_amp[16])
    @test eltype(res.trig_t) == eltype(t)

    # the window starts 5 µs before the trigger, its baseline stretch is at zero and the pulse follows
    n_pre = round(Int, 5.0u"µs" / Δt)
    @test abs(mean(res.pulse[1][1:250])) < 0.1
    @test maximum(res.pulse[1][n_pre:n_pre+40]) ≈ 20.0 atol = 4.0
    @test abs(mean(res.noise[16])) < 0.1

    # the partitioned method: the same selection for one partition, one row per waveform for several
    res_part = dsp_sipm_pulse_selection_compressed(length(wvfs), wvfs, config)
    @test res_part.is_single == res.is_single
    @test res_part.pulse[11] == res.pulse[11]
    @test length(dsp_sipm_pulse_selection_compressed(13, wvfs, config)) == length(wvfs)

    # the pulse condition: no isolated pulse stays above the pulse threshold for 2 µs; without it waveform 27 is one
    config_tot = deepcopy(config)
    config_tot.trigger.min_tot_pulse = 2.0u"µs"
    @test !any(dsp_sipm_pulse_selection_compressed(wvfs, config_tot).is_single)
    config_low = deepcopy(config)
    config_low.trigger.n_σ_pulse = 4.0
    @test dsp_sipm_pulse_selection_compressed(wvfs, config_low).is_single[27]
    config_low.trigger.n_σ_pulse = 3.0
    @test_throws ArgumentError dsp_sipm_pulse_selection_compressed(wvfs, config_low)

    config_odd = deepcopy(config)
    config_odd.trigger.passes = 3
    @test_throws ArgumentError dsp_sipm_pulse_selection_compressed(wvfs, config_odd)
    config_bl = deepcopy(config)
    config_bl.baseline_window = [-6.0u"µs", -2.0u"µs"]
    @test_throws ArgumentError dsp_sipm_pulse_selection_compressed(wvfs, config_bl)
end

@testset "dsp_sipm_template, dsp_sipm_noise_spectrum, dsp_sipm_kernel" begin
    n = 6250
    Δt = 16.0u"ns"
    t = range(0.0u"ns", step = Δt, length = n)
    rng = MersenneTwister(5)
    i_mid = round(Int, 48.0u"µs" / Δt)
    pulse_shape(a, i0 = i_mid) = [i < i0 ? 0.0 : a * (1 - exp(-(i - i0) / 3)) * exp(-(i - i0) / 300) for i in 1:n]
    tail(a) = [a * exp(-i / 1875) for i in 1:n]
    signals = vcat([3750.0 .+ 0.3 .* randn(rng, n) .+ pulse_shape(20 + 40 * rand(rng)) for _ in 1:300],     #   1-300: one pulse
                   [3750.0 .+ 0.3 .* randn(rng, n) for _ in 1:60],                                          # 301-360: noise
                   [3750.0 .+ 0.3 .* randn(rng, n) .+ pulse_shape(40.0, round(Int, 10.0u"µs" / Δt))],       # 361: one pulse before min_pretrigger
                   [3750.0 .+ 0.3 .* randn(rng, n) .+ pulse_shape(40.0) .+ tail(1.5)],                      # 362: one pulse on the slow tail of an earlier one
                   [3750.0 .+ 0.3 .* randn(rng, n) .+ tail(1.0)])                                           # 363: noise on the slow tail of an earlier pulse
    wvfs = ArrayOfRDWaveforms((fill(t, length(signals)), VectorOfSimilarVectors(reduce(hcat, signals))))
    template_config = PropDict(
        :trigger => PropDict(:mwa_window => 160.0u"ns", :passes => 4, :n_σ_threshold => 4.0, :n_σ_valley => 4.0, :min_tot => 80.0u"ns",
                             :n_σ_pulse => 5.0, :min_tot_pulse => 64.0u"ns", :min_threshold => -2.0, :max_threshold => 2.0),
        :min_pretrigger => 25.0u"µs", :window => [-5.0u"µs", 30.0u"µs"], :baseline_window => [-5.0u"µs", -1.0u"µs"],
        :n_σ_baseline => 3.0, :flat_window => 10.0u"µs", :n_σ_flat => 3.0, :n_σ_noise => 8.0, :max_noise_wvfs => 100,
        :tail_fit_window => [10.0u"µs", 20.0u"µs"], :min_pulses => 100, :min_noise_traces => 20)
    kernel_config = PropDict(:target_fwhm => 75.0u"ns", :tail_start_n_sigma => 5.0, :tail_end => 2.0u"µs", :tail_max => 0.01,
                             :lambda_range => [1e-9, 1.0], :n_lambda => 37, :kernel_half_length => 48.0u"µs")
    dsp_sel = dsp_sipm_pulse_selection_compressed(wvfs, template_config)

    # pulses of the template: not the early pulse and not the pulse on a slow tail, both isolated pulses of the selection
    tmpl = dsp_sipm_template(dsp_sel, Δt, template_config)
    @test all(dsp_sel.is_single[[361, 362]]) && !any(tmpl.is_pulse[[361, 362]])
    @test tmpl.is_pulse == tmpl.is_pulse .& dsp_sel.is_single
    @test tmpl.n_pulses == count(tmpl.is_pulse) >= 290

    # template: height 1, the decay of the pulses, their height per trigger amplitude, the decay time of the tail
    @test tmpl.n_pretrigger == 312 && tmpl.n_samples == length(tmpl.pulse) == 2188 && tmpl.dt == Δt
    @test maximum(tmpl.pulse) == 1
    i_peak = argmax(tmpl.pulse)
    @test tmpl.n_pretrigger < i_peak < tmpl.n_pretrigger + 30
    @test tmpl.pulse[i_peak+600] / tmpl.pulse[i_peak+300] ≈ exp(-1) rtol = 0.05
    @test tmpl.peak ≈ median(maximum.(dsp_sel.pulse[tmpl.is_pulse]) ./ dsp_sel.trig_amp[tmpl.is_pulse]) rtol = 0.05
    @test tmpl.tail_tau ≈ 300Δt rtol = 0.05
    @test tmpl.tail_level ≈ tmpl.pulse[end] atol = 0.002

    # short pulses that have decayed in the tail fit window: no tail
    short = [k < 313 ? 0.0 : exp(-(k - 313) / 30) for k in 1:2188]
    sel_short = (is_single = trues(200), pretrigger = fill(48.0u"µs", 200), flatness = 0.02 .* randn(rng, 200), trig_amp = fill(2.0, 200),
                 pulse = [10 .* short .+ 0.1 .* randn(rng, 2188) for _ in 1:200])
    flat = dsp_sipm_template(sel_short, Δt, template_config)
    @test iszero(flat.tail_tau) && iszero(flat.tail_level)
    @test flat.peak ≈ 5 rtol = 0.01
    @test maximum(abs, flat.pulse .- short) < 0.01

    # noise spectrum of white noise of σ = 0.3, without the trace on a slow tail
    noise = dsp_sipm_noise_spectrum(dsp_sel, Δt, template_config)
    @test dsp_sel.is_noise[363] && !isempty(dsp_sel.noise[363]) && !noise.is_trace[363]
    @test noise.is_trace == noise.is_trace .& dsp_sel.is_noise
    @test count(dsp_sel.is_noise) - 5 <= noise.n_traces == count(noise.is_trace) < count(dsp_sel.is_noise)
    @test noise.n_traces >= 50
    @test length(noise.psd) == length(noise.freq) == n ÷ 2 + 1
    @test sum(noise.psd) * step(noise.freq) ≈ 0.09 rtol = 0.05

    # kernel: template continued with the tail, target with sum 1 at the trigger, response on the target in the tail window
    kern = dsp_sipm_kernel(tmpl, noise, kernel_config)
    @test length(kern.template) == length(kern.target) == length(kern.response) == n
    @test kern.template[1:2188] == tmpl.pulse
    @test kern.template[2189] ≈ tmpl.tail_level * exp(-Δt / tmpl.tail_tau)
    kern_flat = dsp_sipm_kernel(flat, noise, kernel_config)
    @test kern_flat.template[1:2188] == flat.pulse && all(iszero, kern_flat.template[2189:end])
    @test sum(kern.target) ≈ 1 && argmax(kern.target) == tmpl.n_pretrigger + 1
    @test kern.n_samples == length(kern.coeffs) == 6001 && kern.dt == Δt && kern.target_fwhm == 75.0u"ns"
    @test first(kern.tail) == 10 && last(kern.tail) == 125
    @test length(kern.scan.λ) == 37 && kern.lambda in kern.scan.λ
    @test kern.tail_deviation <= 0.01 && kern.tail_sum <= 0.01
    @test maximum(abs, kern.response[tmpl.n_pretrigger+1 .+ kern.tail] .- kern.target[tmpl.n_pretrigger+1 .+ kern.tail]) ≈ kern.tail_deviation * maximum(kern.target)
    @test argmax(kern.response) == tmpl.n_pretrigger + 1
    @test sum(kern.response) ≈ 1 rtol = 0.05

    # spectra the kernel is built from: zero frequency holds the sums of template, target and coefficients
    @test kern.spectrum.freq == noise.freq
    @test length(kern.spectrum.template) == length(kern.spectrum.target) == length(kern.spectrum.noise) == length(kern.spectrum.gain) == n ÷ 2 + 1
    @test first(kern.spectrum.template) ≈ sum(kern.template)^2 && first(kern.spectrum.target) ≈ 1
    @test first(kern.spectrum.gain) ≈ abs(sum(kern.coeffs))
    @test median(kern.spectrum.noise[3:end-2]) ≈ n * 0.09 rtol = 0.1

    @test_throws ArgumentError dsp_sipm_template(dsp_sel, Δt, merge(template_config, PropDict(:min_pulses => 1000)))
    @test_throws ArgumentError dsp_sipm_template(dsp_sel, Δt, merge(template_config, PropDict(:flat_window => 25.0u"µs")))
    @test_throws ArgumentError dsp_sipm_template(dsp_sel, Δt, merge(template_config, PropDict(:tail_fit_window => [20.0u"µs", 40.0u"µs"])))
    @test_throws ArgumentError dsp_sipm_noise_spectrum(dsp_sel, Δt, merge(template_config, PropDict(:min_noise_traces => 100)))
    @test_throws ArgumentError dsp_sipm_kernel(merge(tmpl, (dt = 8.0u"ns",)), noise, kernel_config)
end
