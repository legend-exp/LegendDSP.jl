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
