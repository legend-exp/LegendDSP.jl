# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    dsp_sipm_pulse_selection_compressed(wvfs::ArrayOfRDWaveforms, config::PropDict)

Candidates for the template and the noise spectrum of a SiPM channel: waveforms with one isolated pulse and waveforms without any pulse, without pole-zero correction and without pars.

# Arguments
- `wvfs::ArrayOfRDWaveforms`: Array of RDWaveforms
- `config::PropDict`: template configuration with `trigger` (`mwa_window`, `passes` (even), `n_σ_threshold`, `n_σ_valley`, `min_tot`, `n_σ_pulse`, `min_tot_pulse`, `min_threshold`, `max_threshold`), `window`, `baseline_window`, `n_σ_baseline`, `flat_window`, `n_σ_noise` and `max_noise_wvfs`

# Output columns
- `n_trig`: number of positive triggers
- `trig_t`, `trig_amp`: time and amplitude on the trigger trace of the first positive trigger, `NaN` without one
- `threshold`: MAD-σ of the trigger trace
- `pretrigger`: time from the start of the trace to the first positive trigger, `NaN` without one
- `is_single`: the waveform holds one isolated pulse with `window` inside the trace
- `is_noise`: the waveform holds no trigger and no sample farther than `n_σ_noise` σ from its baseline
- `flatness`: how much higher the trace starts than the pulse or the trace ends, in ADC, `NaN` for all other waveforms
- `pulse`: baseline-subtracted `window` around the trigger of the `is_single` waveforms, an empty vector otherwise
- `noise`: the full baseline-subtracted trace of the first `max_noise_wvfs` `is_noise` waveforms, an empty vector otherwise

# Processing
- *Trigger trace:* the derivative of the waveform (`DerivativeFilter(1)`) after `passes` moving averages of `mwa_window` (`MovingWindowMultiFilter`). Baseline and noise width come from [`madstats`](@ref) of the samples between `min_threshold` and `max_threshold`, and the baseline is subtracted.
- *Triggers:* an `IntersectExtrema(min_tot)` with the threshold `n_σ_threshold` and the valley `n_σ_valley` times the noise width. Triggers within `2 passes` times `mwa_window` of either trace end are dropped. These triggers count the pulses of a waveform.
- *Isolated pulse:* exactly one positive trigger with `window` around it inside the trace. Its region must hold a trigger of a second `IntersectExtrema(min_tot_pulse)` with the threshold `n_σ_pulse` times the noise width, so that the one pulse stands clearly above the noise. The window is aligned at the trigger time, rounded to the sample, and its baseline, from [`sigmaclipstats`](@ref) of the samples in `baseline_window` with `n_σ_baseline`, is subtracted.
- *Noise trace:* no trigger of either sign. Baseline and noise width come from `sigmaclipstats` of the whole trace with `n_σ_baseline`.
- *Flatness:* the trigger does not see the slow tail of a pulse before the start of the trace. `flatness` is, for an isolated pulse, the mean of the first `flat_window` of the trace minus the baseline of the pulse, and for a noise trace the mean of its first minus the mean of its last `flat_window`. It has a meaning for a pulse only if the pulse comes after the first `flat_window`; the cuts on `pretrigger` and `flatness` are made where the pulses and the noise traces are used.
"""
function dsp_sipm_pulse_selection_compressed end
export dsp_sipm_pulse_selection_compressed

function dsp_sipm_pulse_selection_compressed(wvfs::ArrayOfRDWaveforms, config::PropDict)
    # get dsp meta parameters
    mwa_window        = config.trigger.mwa_window
    mwa_passes        = config.trigger.passes
    n_σ_threshold     = config.trigger.n_σ_threshold
    n_σ_valley        = config.trigger.n_σ_valley
    min_tot           = config.trigger.min_tot
    n_σ_pulse         = config.trigger.n_σ_pulse
    min_tot_pulse     = config.trigger.min_tot_pulse
    min_threshold     = config.trigger.min_threshold
    max_threshold     = config.trigger.max_threshold
    window            = config.window
    baseline_window   = config.baseline_window
    n_σ_baseline      = config.n_σ_baseline
    flat_window       = config.flat_window
    n_σ_noise         = config.n_σ_noise
    max_noise_wvfs    = config.max_noise_wvfs
    iseven(mwa_passes) || throw(ArgumentError("The trigger trace needs an even number of passes to stay centred, got $mwa_passes"))
    n_σ_pulse >= n_σ_threshold || throw(ArgumentError("The pulse threshold n_σ_pulse = $n_σ_pulse must not lie below the trigger threshold n_σ_threshold = $n_σ_threshold"))
    first(window) <= first(baseline_window) < last(baseline_window) <= last(window) || throw(ArgumentError("The baseline window $baseline_window must lie inside the window $window"))

    # get waveform data
    wvfs = decode_data(wvfs)

    # shift waveform by 0 to get Float64 conversion
    wvfs = shift_waveform.(wvfs, 0.0)

    # trigger trace: derivative plus moving-window shaping, baseline subtracted
    wvfs_deriv = DerivativeFilter(1).(wvfs)
    mwa_filter = MovingWindowMultiFilter(mwa_window; passes = mwa_passes)
    wvfs_mwa = mwa_filter.(wvfs_deriv)
    bl_mwa = madstats.(wvfs_mwa, min_threshold, max_threshold)
    wvfs_mwa = shift_waveform.(wvfs_mwa, [-bl.median for bl in bl_mwa])
    thresholds = [bl.sigma for bl in bl_mwa]

    # trigger finding, without the shaping transients at both ends of the trace
    t = first(wvfs).time
    dt = step(t)
    intflt = IntersectExtrema(min_tot)
    mwa_margin = 2 * mwa_passes * mwa_window
    trigs = map(intflt.(wvfs_mwa, n_σ_threshold .* thresholds, n_σ_valley .* thresholds)) do s
        k = findall(x -> first(t) + mwa_margin <= x <= last(t) - mwa_margin, s.t)
        (t = s.t[k], amp = s.amp[k], t_start = s.t_start[k], t_end = s.t_end[k])
    end
    n_trig = [count(>(0), s.amp) for s in trigs]
    first_trig = [findfirst(>(0), s.amp) for s in trigs]
    trig_t = [isnothing(p) ? oftype(float(first(t)), NaN * oneunit(first(t))) : s.t[p] for (s, p) in zip(trigs, first_trig)]
    trig_amp = [isnothing(p) ? NaN : s.amp[p] for (s, p) in zip(trigs, first_trig)]
    intflt_pulse = IntersectExtrema(min_tot_pulse)
    is_pulse = map(trigs, first_trig, intflt_pulse.(wvfs_mwa, n_σ_pulse .* thresholds, n_σ_valley .* thresholds)) do s, p, s_pulse
        !isnothing(p) && any(x -> s.t_start[p] <= x <= s.t_end[p], s_pulse.t)
    end

    # template window and its baseline stretch in samples
    n_samples = length(t)
    n_pre   = round(Int, uconvert(NoUnits, -first(window) / dt))
    n_post  = round(Int, uconvert(NoUnits, last(window) / dt))
    bl_idx  = round(Int, uconvert(NoUnits, (first(baseline_window) - first(window)) / dt)) + 1 : round(Int, uconvert(NoUnits, (last(baseline_window) - first(window)) / dt)) + 1
    n_flat  = round(Int, uconvert(NoUnits, flat_window / dt))

    # isolated pulses: exactly one positive trigger, above the pulse threshold and with the window inside the trace
    trig_idx = [isnan(x) ? 0 : round(Int, uconvert(NoUnits, (x - first(t)) / dt)) + 1 for x in trig_t]
    is_single = [n == 1 && p && i - n_pre >= 1 && i + n_post <= n_samples for (n, p, i) in zip(n_trig, is_pulse, trig_idx)]

    # baseline-subtracted template window (isolated pulses), the first max_noise_wvfs noise traces (no trigger and no sample
    # beyond n_σ_noise σ) and the flatness of both: the start of the trace against the baseline of the pulse or against its end
    pulse = Vector{Vector{Float64}}(undef, length(wvfs))
    noise = Vector{Vector{Float64}}(undef, length(wvfs))
    is_noise = falses(length(wvfs))
    flatness = fill(NaN, length(wvfs))
    n_noise = 0
    for k in eachindex(wvfs)
        y = wvfs[k].signal
        pulse[k] = Float64[]
        noise[k] = Float64[]
        if is_single[k]
            i = trig_idx[k]
            y_window = view(y, i-n_pre:i+n_post)
            bl = sigmaclipstats(view(y_window, bl_idx), n_σ_baseline).mean
            pulse[k] = y_window .- bl
            flatness[k] = mean(view(y, 1:n_flat)) - bl
        elseif isempty(trigs[k].t)
            bl = sigmaclipstats(y, n_σ_baseline)
            is_noise[k] = maximum(x -> abs(x - bl.mean), y) < n_σ_noise * bl.sigma
            is_noise[k] || continue
            flatness[k] = mean(view(y, 1:n_flat)) - mean(view(y, n_samples-n_flat+1:n_samples))
            n_noise < max_noise_wvfs && (noise[k] = y .- bl.mean)
            n_noise += 1
        end
    end

    pretrigger = trig_t .- first(t)

    # output Table
    return TypedTables.Table(
        n_trig = n_trig,
        trig_t = trig_t,
        trig_amp = trig_amp,
        threshold = thresholds,
        pretrigger = pretrigger,
        is_single = is_single,
        is_noise = is_noise,
        flatness = flatness,
        pulse = VectorOfVectors(pulse),
        noise = VectorOfVectors(noise)
    )
end

function dsp_sipm_pulse_selection_compressed(n_max_wvfs::Int, wvfs::ArrayOfRDWaveforms, config::PropDict)
    # flatten over waveform partitions, the noise trace budget shared between them
    n_parts = cld(length(wvfs), n_max_wvfs)
    config_part = merge(config, PropDict(:max_noise_wvfs => cld(config.max_noise_wvfs, n_parts)))
    LegendDataTypes.fast_flatten([dsp_sipm_pulse_selection_compressed(wf, config_part) for wf in Iterators.partition(wvfs, n_max_wvfs)])
end
