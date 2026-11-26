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
- *Flatness:* the trigger does not see the slow tail of a pulse before the start of the trace. `flatness` is, for an isolated pulse, the mean of the first `flat_window` of the trace minus the baseline of the pulse, and for a noise trace the mean of its first minus the mean of its last `flat_window`. It has a meaning for a pulse only if the pulse comes after the first `flat_window`; the cuts on `pretrigger` and `flatness` are made where the pulses and the noise traces are used ([`dsp_sipm_template`](@ref), [`dsp_sipm_noise_spectrum`](@ref)).
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

"""
    dsp_sipm_template(dsp_sel, dt::Unitful.Time, config::PropDict)

Template of a SiPM channel: the pulse shape of the isolated pulses that [`dsp_sipm_pulse_selection_compressed`](@ref) selected.

# Arguments
- `dsp_sel`: table of `dsp_sipm_pulse_selection_compressed` with the columns `is_single`, `pretrigger`, `flatness`, `trig_amp` and `pulse`
- `dt`: sampling interval of the waveforms
- `config::PropDict`: template configuration with `window`, `baseline_window` and `flat_window` (as used for the selection), `min_pretrigger`, `n_σ_flat`, `tail_fit_window` and `min_pulses`

# Returns
- `pulse`: the template, height 1, on the samples of `window`
- `n_samples`: number of samples of `pulse`
- `n_pretrigger`: number of samples of `pulse` before the trigger
- `dt`: sampling interval
- `peak`: pulse height in ADC per unit of trigger amplitude
- `n_pulses`: number of pulses the template is built from
- `is_pulse`: which rows of `dsp_sel` these are
- `tail_tau`, `tail_level`: decay time of the tail of the template and its level at the end of `pulse`, both zero if the tail does not decay in `tail_fit_window`

All of these are fields of one named tuple.

# Processing
- *Pulses:* the isolated pulses that come at least `min_pretrigger` after the start of the trace, so that the first `flat_window` of the trace ends before their baseline window, and whose `flatness` lies within `n_σ_flat` MAD-σ of the median flatness of these pulses ([`madstats`](@ref)): a pulse on the slow tail of an earlier one is left out.
- *Template:* [`wiener_template`](@ref) of these pulses with their trigger amplitudes: every pulse divided by its trigger amplitude, the median per sample, divided by its maximum. Throws below `min_pulses` pulses.
- *Tail:* [`tailstats`](@ref) of the template in `tail_fit_window`, given relative to the trigger, for a template that has not decayed at the end of `window`.
"""
function dsp_sipm_template end
export dsp_sipm_template

function dsp_sipm_template(dsp_sel, dt::Unitful.Time{<:Real}, config::PropDict)
    # get meta parameters
    window              = config.window
    baseline_window     = config.baseline_window
    flat_window         = config.flat_window
    min_pretrigger      = config.min_pretrigger
    n_σ_flat            = config.n_σ_flat
    tail_fit_window     = config.tail_fit_window
    min_pulses          = config.min_pulses
    min_pretrigger + first(baseline_window) >= flat_window || throw(ArgumentError("The flat window $flat_window must end before the baseline window $baseline_window of a pulse at min_pretrigger = $min_pretrigger"))
    first(window) <= first(tail_fit_window) < last(tail_fit_window) <= last(window) || throw(ArgumentError("The tail fit window $tail_fit_window must lie inside the window $window"))

    # isolated pulses late enough in the trace and on a flat trace, one column each (QC)
    is_pulse = dsp_sel.is_single .& (dsp_sel.pretrigger .>= min_pretrigger)
    flat = madstats(dsp_sel.flatness[is_pulse])
    is_pulse = is_pulse .& (abs.(dsp_sel.flatness .- flat.median) .<= n_σ_flat * flat.sigma)
    n_pulses = count(is_pulse)
    n_pulses >= min_pulses || throw(ArgumentError("Only $n_pulses isolated pulses, at least $min_pulses are needed for the template"))
    windows = stack(dsp_sel.pulse[is_pulse])
    n_pretrigger = round(Int, uconvert(NoUnits, -first(window) / dt))

    # template: median of the pulses over their trigger amplitude, height 1
    (; template, peak) = wiener_template(windows, dsp_sel.trig_amp[is_pulse])
    n_samples = length(template)

    # tail: decay time and level at the end of the template
    i_from, i_to = n_pretrigger + 1 .+ round.(Int, uconvert.(NoUnits, (first(tail_fit_window), last(tail_fit_window)) ./ dt))
    tail = tailstats(template, i_from, i_to)
    has_tail = tail.τ > zero(tail.τ)
    tail_tau = has_tail ? tail.τ * dt : zero(tail.τ * dt)
    tail_level = has_tail ? exp(tail.mean - (n_samples - (i_from + i_to) / 2) / tail.τ) : zero(eltype(template))

    return (
        pulse = template,
        n_samples = n_samples,
        n_pretrigger = n_pretrigger,
        dt = dt,
        peak = peak,
        n_pulses = n_pulses,
        is_pulse = is_pulse,
        tail_tau = tail_tau,
        tail_level = tail_level
    )
end

"""
    dsp_sipm_noise_spectrum(dsp_sel, dt::Unitful.Time, config::PropDict)

Noise spectrum of a SiPM channel from the noise traces that [`dsp_sipm_pulse_selection_compressed`](@ref) kept.

# Arguments
- `dsp_sel`: table of `dsp_sipm_pulse_selection_compressed` with the columns `is_noise`, `flatness` and `noise`
- `dt`: sampling interval of the waveforms
- `config::PropDict`: template configuration with `n_σ_flat` and `min_noise_traces`

# Returns
- `freq`, `psd`: frequencies and one-sided noise power per frequency, in ADC² times the unit of `dt`
- `n_traces`: number of noise traces the spectrum is built from
- `is_trace`: which rows of `dsp_sel` these are

All of these are fields of one named tuple.

# Processing
- *Noise traces:* the kept traces whose `flatness` lies within `n_σ_flat` MAD-σ of the median flatness of all `is_noise` waveforms ([`madstats`](@ref)): a trace on the slow tail of an earlier pulse is left out.
- *Noise spectrum:* [`wiener_noise_spectrum`](@ref) of these traces. Throws below `min_noise_traces` traces.
"""
function dsp_sipm_noise_spectrum end
export dsp_sipm_noise_spectrum

function dsp_sipm_noise_spectrum(dsp_sel, dt::Unitful.Time{<:Real}, config::PropDict)
    # get meta parameters
    n_σ_flat            = config.n_σ_flat
    min_noise_traces    = config.min_noise_traces

    # kept noise traces on a flat trace, one column each (QC)
    flat = madstats(dsp_sel.flatness[dsp_sel.is_noise])
    is_trace = .!isempty.(dsp_sel.noise) .& (abs.(dsp_sel.flatness .- flat.median) .<= n_σ_flat * flat.sigma)
    n_traces = count(is_trace)
    n_traces >= min_noise_traces || throw(ArgumentError("Only $n_traces noise traces, at least $min_noise_traces are needed for the noise spectrum"))
    traces = stack(dsp_sel.noise[is_trace])

    noise = wiener_noise_spectrum(traces, dt)

    return (
        freq = noise.freq,
        psd = noise.psd,
        n_traces = n_traces,
        is_trace = is_trace
    )
end

"""
    dsp_sipm_kernel(template, noise, config::PropDict)

Deconvolution kernel of a SiPM channel: turns every pulse of the shape of `template` into a narrow Gaussian whose area is the pulse height in ADC.

# Arguments
- `template`: result of [`dsp_sipm_template`](@ref)
- `noise`: result of [`dsp_sipm_noise_spectrum`](@ref)
- `config::PropDict`: kernel configuration with `target_fwhm`, `tail_start_n_sigma`, `tail_end`, `tail_max`, `lambda_range`, `n_lambda` and `kernel_half_length`

# Returns
- `coeffs`: the kernel, lag 0 at its middle sample; `n_samples`: its number of samples, odd; `dt`: sampling interval
- `lambda`: the chosen noise weight; `tail_deviation`, `tail_sum`: largest deviation and largest running sum of the response from the target in the tail window, relative to the maximum and the sum of the target
- `noise_out`: standard deviation of a sample of the filtered noise in ADC
- `target_fwhm`
- `template`, `target`, `response`: the continued template, the target and the template through the kernel, on one axis as long as the noise traces with `template.n_pretrigger` samples before the trigger
- `tail`: lags after the trigger, in samples, at which response and target are compared
- `scan`: `λ`, `tail_deviation`, `tail_sum` and `noise_out` for every noise weight of the scan
- `spectrum`: on the frequencies `freq` of the noise spectrum, the squared magnitudes of the Fourier transforms of the continued template (`template`) and the target (`target`), the expected one of a noise trace (`noise`) and the `gain` of the kernel

All of these are fields of one named tuple.

# Processing
- *Template:* the template, continued to the length of the noise traces with the exponential of `template.tail_tau` and `template.tail_level`, or with zeros if `tail_tau` is zero. Throws if template and noise spectrum differ in their sampling interval.
- *Target:* a Gaussian of `target_fwhm` and sum 1 at the trigger.
- *Kernel:* [`wiener_optimization`](@ref) with `n_lambda` logarithmic steps of the noise weight over `lambda_range`, the tail window from `tail_start_n_sigma` σ of the target to `tail_end` after the trigger with the limit `tail_max`, and `kernel_half_length` on either side.
"""
function dsp_sipm_kernel end
export dsp_sipm_kernel

function dsp_sipm_kernel(template, noise, config::PropDict)
    # get meta parameters
    target_fwhm         = config.target_fwhm
    tail_start_n_sigma  = config.tail_start_n_sigma
    tail_end            = config.tail_end
    tail_max            = config.tail_max
    lambda_range        = config.lambda_range
    n_lambda            = config.n_lambda
    kernel_half_length  = config.kernel_half_length

    # grid of the noise traces
    dt = template.dt
    n_pre = template.n_pretrigger
    m, n = length(template.pulse), 2 * (length(noise.psd) - 1)
    inv(n * step(noise.freq)) ≈ dt || throw(ArgumentError("The noise spectrum is sampled at $(inv(n * step(noise.freq))), the template at $dt"))

    # the template continued with its tail to the length of the noise traces
    T = eltype(template.pulse)
    continuation = iszero(template.tail_tau) ? zeros(T, n - m) : template.tail_level .* exp.(-(1:n-m) .* uconvert(NoUnits, dt / template.tail_tau))
    template_long = vcat(template.pulse, continuation)

    # target: Gaussian with sum 1 at the trigger
    σ_target = uconvert(NoUnits, target_fwhm / dt) / (2 * sqrt(2 * log(2)))
    target = exp.(-((1:n) .- (n_pre + 1)) .^ 2 ./ (2 * σ_target^2))
    target ./= sum(target)

    # kernel with the largest noise weight that keeps the tail
    lambdas = exp10.(range(log10(first(lambda_range)), log10(last(lambda_range)), length = n_lambda))
    tail = ceil(Int, tail_start_n_sigma * σ_target):round(Int, uconvert(NoUnits, tail_end / dt))
    n_half = round(Int, uconvert(NoUnits, kernel_half_length / dt))
    kernel = wiener_optimization(noise, template_long, target, n_pre; lambdas, tail, tail_max, n_half)
    i_λ = findfirst(==(kernel.λ), lambdas)
    n_samples = length(kernel.coeffs)
    tail_deviation = kernel.scan.tail_deviation[i_λ]
    tail_sum = kernel.scan.tail_sum[i_λ]

    return (
        coeffs = kernel.coeffs,
        n_samples = n_samples,
        dt = dt,
        lambda = kernel.λ,
        tail_deviation = tail_deviation,
        tail_sum = tail_sum,
        noise_out = kernel.noise_out,
        target_fwhm = target_fwhm,
        template = template_long,
        target = target,
        response = kernel.response,
        tail = tail,
        scan = kernel.scan,
        spectrum = kernel.spectrum
    )
end
