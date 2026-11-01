# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    dsp_sipm_template_compressed(wvfs::ArrayOfRDWaveforms, config::PropDict)

Waveform selection for the SiPM template extraction: derivative, moving-window shaping and a MAD threshold
per waveform, without pole-zero correction and without pars. A waveform `is_single` if it has exactly one
trigger, the template window lies inside the trace and at least as many baseline samples as the window
has before the trigger precede it; the window is aligned at the threshold crossing of the shaped trace
(rounded to the sample) and returned baseline subtracted as `pulse`. A
waveform `is_noise` if it has no trigger and no sample farther than `n_σ_noise` MAD-σ from its baseline (the
trigger-free condition first: the MAD-σ of a trace full of pulses is inflated by their tails); the first
`max_noise_wvfs` of them are returned baseline subtracted in full as `noise` (the others as empty vectors).
The baseline is the mean of the samples (before the window, or of the whole trace) within `n_σ_baseline` MAD-σ of
their median.

# Arguments
- `wvfs::ArrayOfRDWaveforms`: Array of RDWaveforms
- `config::PropDict`: template configuration (`trigger`, `window`, `n_σ_baseline`, `n_σ_noise`, `max_noise_wvfs`)

# Returns
- `Table`: `n_trig`, `trig_pos`, `trig_max`, `threshold`, `is_single`, `is_noise`, `pulse`, `noise`
"""
function dsp_sipm_template_compressed end
export dsp_sipm_template_compressed

function dsp_sipm_template_compressed(wvfs::ArrayOfRDWaveforms, config::PropDict)
    # get dsp meta parameters
    mwa_window        = config.trigger.mwa_window
    n_σ_threshold     = config.trigger.n_σ_threshold
    min_threshold     = config.trigger.min_threshold
    max_threshold     = config.trigger.max_threshold
    min_tot_intersect = config.trigger.min_tot_intersect
    max_tot_intersect = config.trigger.max_tot_intersect
    window            = config.window
    n_σ_baseline      = config.n_σ_baseline
    n_σ_noise         = config.n_σ_noise
    max_noise_wvfs    = config.max_noise_wvfs

    # get waveform data
    wvfs = decode_data(wvfs)

    # shift waveform by 0 to get Float64 conversion
    wvfs = shift_waveform.(wvfs, 0.0)

    # derivative plus moving-window shaping: pulses become bumps, tails vanish
    wvfs_deriv = DerivativeFilter(1).(wvfs)
    mwa_filter = MovingWindowMultiFilter(mwa_window)
    wvfs_mwa = mwa_filter.(wvfs_deriv)

    # trigger finding with a MAD threshold per waveform
    intflt = IntersectMaximum(min_tot_intersect, max_tot_intersect)
    thresholds = thresholdstats_mad.(wvfs_mwa, min_threshold, max_threshold)
    inters = intflt.(wvfs_mwa, n_σ_threshold .* thresholds)
    n_trig = length.(inters.x)

    # template window in samples
    t = first(wvfs).time
    dt = step(t)
    n_samples = length(t)
    n_pre   = round(Int, uconvert(NoUnits, -first(window) / dt))
    n_post  = round(Int, uconvert(NoUnits, last(window) / dt))

    # single-trigger waveforms, aligned at the threshold crossing of the shaped trace
    trig_idx = [isempty(x) ? 0 : round(Int, uconvert(NoUnits, (first(x) - first(t)) / dt)) + 1 for x in inters.x]
    is_single = [n == 1 && i - 2n_pre >= 1 && i + n_post <= n_samples for (n, i) in zip(n_trig, trig_idx)]

    # baseline: mean of the samples within n_σ_baseline MAD-σ of their median, and that σ
    function baseline(seg)
        med = median(seg)
        σ = 1.4826 * median(abs.(seg .- med))
        (mean(v for v in seg if abs(v - med) <= n_σ_baseline * σ), σ)
    end

    # baseline-subtracted template window (single-trigger waveforms) and full trace (trigger-free waveforms with no
    # sample beyond n_σ_noise MAD-σ, the first max_noise_wvfs of them)
    pulse = Vector{Vector{Float64}}(undef, length(wvfs))
    noise = Vector{Vector{Float64}}(undef, length(wvfs))
    is_noise = falses(length(wvfs))
    n_noise = 0
    for k in eachindex(wvfs)
        y = wvfs[k].signal
        if is_single[k]
            i = trig_idx[k]
            pulse[k] = y[i-n_pre:i+n_post] .- first(baseline(view(y, 1:i-n_pre-1)))
        else
            pulse[k] = Float64[]
        end
        noise[k] = Float64[]
        if n_trig[k] == 0
            bl, σ = baseline(y)
            y_bl = y .- bl
            is_noise[k] = maximum(abs, y_bl) < n_σ_noise * σ
            if is_noise[k] && n_noise < max_noise_wvfs
                n_noise += 1
                noise[k] = y_bl
            end
        end
    end

    # output Table
    return TypedTables.Table(
        n_trig = n_trig,
        trig_pos = [isempty(x) ? NaN * oneunit(first(t)) : first(x) for x in inters.x],
        trig_max = [isempty(m) ? NaN : first(m) for m in inters.max],
        threshold = thresholds,
        is_single = is_single,
        is_noise = is_noise,
        pulse = VectorOfVectors(pulse),
        noise = VectorOfVectors(noise)
    )
end

function dsp_sipm_template_compressed(n_max_wvfs::Int, wvfs::ArrayOfRDWaveforms, config::PropDict)
    # flatten over waveform partitions, the noise trace budget shared between them
    n_parts = cld(length(wvfs), n_max_wvfs)
    config_part = merge(config, PropDict(:max_noise_wvfs => cld(config.max_noise_wvfs, n_parts)))
    LegendDataTypes.fast_flatten([dsp_sipm_template_compressed(wf, config_part) for wf in Iterators.partition(wvfs, n_max_wvfs)])
end
