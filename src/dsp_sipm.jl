# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    dsp_sipm(data::Table, config::PropDict, pars::PropDict)

DSP chain for SiPM waveforms. Returns a table with one row per waveform.

# Input columns
`waveform`, `baseline`, `timestamp`, `eventnumber`, `daqenergy`

# Arguments
- `data`: table of waveforms; encoded columns are decoded with `decode_data`
- `config`: `PropDict` with the entries `sg_flt_degree`, `t0_hpge_window` (start and stop of the window in which the HPGe signal is expected, used for the `_lar` parameters) and `filters`. `filters.sg` holds the settings of the Savitzky-Golay pipeline and `filters.trap` those of the trapezoidal pipeline. Both contain `min_tot_intersect`, `max_tot_intersect`, `min_threshold`, `max_threshold`, `n_σ_threshold`, `min_dc_threshold`, `max_dc_threshold` and `n_σ_dc_threshold`; `filters.trap` also contains `rt`, `ft` and `pz_tau`. `filters.xt_tag` holds the settings of the cross-talk tag pipeline: `window_length`, `n_σ_threshold`, `min_tot_intersect` and `max_tot_intersect`. `filters.mwa` holds the settings of the MWA pipeline: `wl_trigger`, `passes` (even), `n_σ_threshold`, `min_threshold`, `max_threshold`, `solve_overlap` and `merge_distance`.
- `pars`: `PropDict` with the optimized Savitzky-Golay window length `pars.sg.wl`, the MWA energy shaping length `pars.mwa.wl` and, when the template step has run for the channel, its inverse kernel `pars.kernel` (`taps`, `dt`)

# Output columns
Times are in the unit of the waveform time axis unless noted.
- `blfc`, `timestamp`, `eventID_fadc`, `e_fc`: baseline, timestamp, event number and energy from the FADC
- `t_max`, `t_min`: times of the waveform maximum and minimum, in µs
- `t_max_lar`, `t_min_lar`: the same within `t0_hpge_window`, in µs
- `e_max`, `e_min`: waveform maximum and minimum
- `e_max_lar`, `e_min_lar`: the same within `t0_hpge_window`
- `blmean`, `blsigma`, `blslope`, `bloffset`: statistics of the integrated Savitzky-Golay derivative from the start of the waveform to the first Savitzky-Golay trigger (to the end of the waveform if there is none, but at least three samples)
- `wfmean`, `wfsigma`, `wfslope`, `wfoffset`: the same statistics over the whole integrated waveform
- `threshold_sg`: trigger threshold of the Savitzky-Golay pipeline
- `threshold_DC_sg`: discharge threshold of the Savitzky-Golay pipeline
- `trig_pos_sg`, `trig_max_sg`: times and maxima of the triggers, one vector per waveform
- `trig_pos_DC_sg`, `trig_max_DC_sg`: times and maxima of the discharge triggers
- `threshold_trap`, `threshold_DC_trap`: trigger and discharge thresholds of the trapezoidal pipeline
- `trig_pos_trap`, `trig_pos_high_trap`, `trig_pos_tot_trap`, `trig_max_trap`: time, falling-edge time, time over threshold and maximum of the trapezoidal triggers
- `trig_pos_DC_trap`, `trig_pos_high_DC_trap`, `trig_pos_tot_DC_trap`, `trig_max_DC_trap`: the same for the trapezoidal discharge triggers
- `threshold_xt_tag`: threshold of the cross-talk tag pipeline
- `trig_pos_xt_tag`, `trig_pos_high_xt_tag`, `trig_pos_tot_xt_tag`, `trig_max_xt_tag`: start time, end time, time over threshold and maximum of each region of the long-window integrated waveform that is above the threshold
- `threshold_mwa`: MAD-σ of the site trace of the MWA pipeline (the trigger trace with `solve_overlap`, the energy trace without)
- `trig_pos_mwa`, `trig_pos_tot_mwa`, `trig_max_mwa`: site positions (every local maximum above `n_σ_threshold` σ and minimum below its negative), time over threshold of the excursion a site lies in, and the signed site amplitudes on the energy shaping, solved for the overlap with `solve_overlap`

# Processing
- The waveform is converted to `Float64`.
- *Savitzky-Golay pipeline:* the waveform is filtered with a `SavitzkyGolayFilter(pars.sg.wl, config.sg_flt_degree, 1)`, which gives the smoothed derivative. The threshold is `n_σ_threshold` times [`thresholdstats_mad`](@ref) of the samples between `min_threshold` and `max_threshold`, and the triggers are found with [`IntersectMaximum`](@ref).
- The derivative is integrated with `IntegratorFilter(gain=1)`, which gives a smoothed waveform, and the baseline statistics are determined on it.
- *Discharge detection:* the median of the integrated waveform is taken as its baseline and subtracted, and the result is multiplied by -1. The thresholds are computed in the same way as above from the `dc` settings of `filters.sg`, and the triggers are found with the Savitzky-Golay `IntersectMaximum`.
- *Trapezoidal pipeline:* the integrated waveform is pole-zero corrected with `pz_tau` and filtered with a `TrapezoidalChargeFilter(rt, ft)`. Thresholds and triggers are computed in the same way as for the Savitzky-Golay pipeline, with the settings of `filters.trap`.
- *Trapezoidal discharge detection:* the thresholds are computed from the flipped, baseline-subtracted integrated waveform with the `dc` settings of `filters.trap`, and the triggers are found with the `IntersectMaximum` of the trapezoidal settings (`filters.trap.min_tot_intersect` and `max_tot_intersect`).
- *Cross-talk tag pipeline:* the waveform is filtered with a `SavitzkyGolayFilter(filters.xt_tag.window_length, config.sg_flt_degree, 1)` and integrated, and the median of the result is subtracted. The threshold is `filters.xt_tag.n_σ_threshold` times [`thresholdstats_mad`](@ref) of all samples, without lower and upper bounds. The regions above the threshold are found with an `IntersectMaximum(filters.xt_tag.min_tot_intersect, filters.xt_tag.max_tot_intersect)`; the waveform is not flipped.
- *MWA pipeline:* the baseline-subtracted (median) waveform goes through the inverse kernel `pars.kernel.taps` of the template step with a `CenteredConvolutionFilter`, or through the derivative kernel `[1, -1]` when the channel has no kernel. The result is shaped with `passes` moving averages of the trigger length `wl_trigger` and of the energy length `pars.mwa.wl` (`MovingWindowMultiFilter`), each scaled so that a unit stick has peak 1. The sites are found with [`extrema_over_threshold`](@ref) on the trigger trace (`solve_overlap`) or on the energy trace, with the threshold `n_σ_threshold` times [`thresholdstats_mad`](@ref) of the samples between `min_threshold` and `max_threshold`. With `solve_overlap` the amplitudes of all sites are solved on the energy trace with the [`SiteSolver`](@ref), sites closer than `merge_distance` times the energy length merging; sites within `2 passes` times the longer shaping length of either trace end are dropped.
"""
function dsp_sipm(data::Q, config::PropDict, pars::PropDict) where {Q <: Table}
    # common config
    sg_flt_degree         = config.sg_flt_degree
    t0_hpge_window        = config.t0_hpge_window

    # get SG config parameters
    sg_config = config.filters.sg
    min_tot_intersect     = sg_config.min_tot_intersect
    max_tot_intersect     = sg_config.max_tot_intersect
    min_threshold         = sg_config.min_threshold
    max_threshold         = sg_config.max_threshold
    n_σ_threshold         = sg_config.n_σ_threshold
    min_dc_threshold      = sg_config.min_dc_threshold
    max_dc_threshold      = sg_config.max_dc_threshold
    n_σ_dc_threshold      = sg_config.n_σ_dc_threshold

    # get trap config parameters
    trap_config = config.filters.trap
    trap_rt                = trap_config.rt
    trap_ft                = trap_config.ft
    trap_pz_tau            = trap_config.pz_tau
    trap_min_tot_intersect = trap_config.min_tot_intersect
    trap_max_tot_intersect = trap_config.max_tot_intersect
    trap_min_threshold     = trap_config.min_threshold
    trap_max_threshold     = trap_config.max_threshold
    trap_n_σ_threshold     = trap_config.n_σ_threshold
    trap_min_dc_threshold  = trap_config.min_dc_threshold
    trap_max_dc_threshold  = trap_config.max_dc_threshold
    trap_n_σ_dc_threshold  = trap_config.n_σ_dc_threshold

    xt_tag_config = config.filters.xt_tag
    xt_tag_window_length     = xt_tag_config.window_length
    xt_tag_n_σ_threshold     = xt_tag_config.n_σ_threshold
    xt_tag_min_tot_intersect = xt_tag_config.min_tot_intersect
    xt_tag_max_tot_intersect = xt_tag_config.max_tot_intersect

    # get MWA config parameters
    mwa_config = config.filters.mwa
    mwa_wl_trigger         = mwa_config.wl_trigger
    mwa_passes             = mwa_config.passes
    mwa_n_σ_threshold      = mwa_config.n_σ_threshold
    mwa_min_threshold      = mwa_config.min_threshold
    mwa_max_threshold      = mwa_config.max_threshold
    mwa_solve_overlap      = mwa_config.solve_overlap
    mwa_merge_distance     = mwa_config.merge_distance
    iseven(mwa_passes) || throw(ArgumentError("The MWA chain needs an even number of passes to stay centred, got $mwa_passes"))

    # unpack optimization parameters and, if the template step has run, the inverse kernel
    sg_window_length = pars.sg.wl
    mwa_wl_energy    = pars.mwa.wl
    has_kernel       = haskey(pars, :kernel)

    # get waveform data
    data = decode_data(data)
    wvfs = data.waveform
    blfc = data.baseline
    ts   = data.timestamp
    evID = data.eventnumber
    efc  = data.daqenergy

    # shift waveform by 0 to get Float64 conversion
    wvfs = shift_waveform.(wvfs, 0.0)

    # get wvf maximum and minimum with timepoints
    estats = extremestats.(wvfs)

    # get wvf maximum and minimum with timepoints for truncated waveform
    uflt_trunc = TruncateFilter(first(t0_hpge_window)..last(t0_hpge_window))
    estats_trunc = extremestats.(uflt_trunc.(wvfs))

    # === SG pipeline ===
    # savitzky golay filter: takes derivative of waveform plus smoothing
    sgflt_savitz = SavitzkyGolayFilter(sg_window_length, sg_flt_degree, 1)
    wvfs_sg = sgflt_savitz.(wvfs)

    # trigger finding on SG-filtered waveforms
    intflt_sg = IntersectMaximum(min_tot_intersect, max_tot_intersect)
    inters_thres = thresholdstats_mad.(wvfs_sg, min_threshold, max_threshold)
    inters = intflt_sg.(wvfs_sg, n_σ_threshold .* inters_thres)

    # integrate derivative
    integrator_filter = IntegratorFilter(gain=1)
    wvfs_int = integrator_filter.(wvfs_sg)

    # get blstats on integrated waveforms
    time_min = minimum(wvfs_int[1].time)
    time_max = last(wvfs_int[1].time)
    # baseline window ends at the first trigger (end of waveform if none), but spans at least 3 samples
    t_bl_min = time_min + 3*step(wvfs_int[1].time)
    t_bl_stop = map(x -> oftype(time_max, max(minimum(x; init = time_max), t_bl_min)), inters.x)
    bl_stats = signalstats.(wvfs_int, time_min, t_bl_stop)
    sigstats = signalstats.(wvfs_int, time_min, time_max)
    
    # flip around x-axis for discharge detection, relative to the baseline (median) of each integrated waveform
    wvfs_int_flip = multiply_waveform.(shift_waveform.(wvfs_int, .-median.(wvfs_int.signal)), -1.0)
    inters_thres_DC = thresholdstats_mad.(wvfs_int_flip, min_dc_threshold, max_dc_threshold)
    inters_DC = intflt_sg.(wvfs_int_flip, n_σ_dc_threshold .* inters_thres_DC)

    # === Trap pipeline ===
    # PZ correction on integrated waveforms
    pz_filter = InvCRFilter(trap_pz_tau)
    wvfs_pz = pz_filter.(wvfs_int)

    # trapezoidal charge filter
    trap_filter = TrapezoidalChargeFilter(trap_rt, trap_ft)
    wvfs_trap = trap_filter.(wvfs_pz)

    # trigger finding on trap-filtered waveforms
    intflt_trap = IntersectMaximum(trap_min_tot_intersect, trap_max_tot_intersect)
    inters_thres_trap = thresholdstats_mad.(wvfs_trap, trap_min_threshold, trap_max_threshold)
    inters_trap = intflt_trap.(wvfs_trap, trap_n_σ_threshold .* inters_thres_trap)

    # discharge detection on flipped integrated waveforms 
    inters_thres_DC_trap = thresholdstats_mad.(wvfs_int_flip, trap_min_dc_threshold, trap_max_dc_threshold)
    inters_DC_trap = intflt_trap.(wvfs_int_flip, trap_n_σ_dc_threshold .* inters_thres_DC_trap)

    wvfs_xt_tag = integrator_filter.(SavitzkyGolayFilter(xt_tag_window_length, sg_flt_degree, 1).(wvfs))
    wvfs_xt_tag = shift_waveform.(wvfs_xt_tag, .-median.(wvfs_xt_tag.signal))
    intflt_xt_tag = IntersectMaximum(xt_tag_min_tot_intersect, xt_tag_max_tot_intersect)
    inters_thres_xt_tag = thresholdstats_mad.(wvfs_xt_tag)
    inters_xt_tag = intflt_xt_tag.(wvfs_xt_tag, xt_tag_n_σ_threshold .* inters_thres_xt_tag)

    # === MWA pipeline ===
    # photon trace: the baseline-subtracted waveform through the inverse kernel of the template step; a channel without
    # template pars takes the derivative (kernel [1, -1]) instead
    t_axis = first(wvfs).time
    dt = step(t_axis)
    !has_kernel || pars.kernel.dt ≈ dt || throw(ArgumentError("The kernel taps are sampled at $(pars.kernel.dt), the waveforms at $dt"))
    photon_filter = CenteredConvolutionFilter(has_kernel ? Float64.(pars.kernel.taps) : [0.0, 1.0, -1.0])
    wvfs_bl = shift_waveform.(wvfs, [-median(wf.signal) for wf in wvfs])
    wvfs_photons = photon_filter.(wvfs_bl)

    # short (trigger) and long (energy) shaping, each scaled so that a unit stick (one photon after the kernel) has peak 1
    mwa_trigger_filter = MovingWindowMultiFilter(mwa_wl_trigger; passes = mwa_passes)
    mwa_energy_filter  = MovingWindowMultiFilter(mwa_wl_energy; passes = mwa_passes)
    stick = RDWaveform(t_axis, [i == length(t_axis) ÷ 2 ? 1.0 : 0.0 for i in eachindex(t_axis)])
    response_trigger = mwa_trigger_filter(stick).signal
    response_energy  = mwa_energy_filter(stick).signal
    wvfs_trigger = multiply_waveform.(mwa_trigger_filter.(wvfs_photons), 1 / maximum(response_trigger))
    wvfs_energy  = multiply_waveform.(mwa_energy_filter.(wvfs_photons), 1 / maximum(response_energy))
    wvfs_trigger = shift_waveform.(wvfs_trigger, [-median(wf.signal) for wf in wvfs_trigger])
    wvfs_energy  = shift_waveform.(wvfs_energy, [-median(wf.signal) for wf in wvfs_energy])

    # sites: the extrema beyond n σ of the trigger trace (matrix method) or of the energy trace
    wvfs_sites = mwa_solve_overlap ? wvfs_trigger : wvfs_energy
    mwa_thresholds = thresholdstats_mad.(wvfs_sites, mwa_min_threshold, mwa_max_threshold)
    sites = extrema_over_threshold.(wvfs_sites, mwa_n_σ_threshold .* mwa_thresholds)

    # matrix method: the amplitudes of all sites solved on the long shaping (otherwise the extrema of the energy trace)
    if mwa_solve_overlap
        solver = SiteSolver(response_energy ./ maximum(response_energy), argmax(response_energy), mwa_merge_distance * uconvert(NoUnits, mwa_wl_energy / dt))
        sites = map(zip(wvfs_energy, sites)) do (wf, s)
            sol = solver(wf.signal, uconvert.(NoUnits, (s.x .- first(t_axis)) ./ dt) .+ 1, uconvert.(NoUnits, s.x_tot ./ dt), Int.(sign.(s.max)))
            (x = convert.(eltype(t_axis), first(t_axis) .+ (sol.pos .- 1) .* dt), x_tot = convert.(eltype(t_axis), sol.tot .* dt), max = sol.A)
        end
    end

    # without the shaping transients at both ends of the trace
    mwa_margin = 2 * mwa_passes * max(mwa_wl_trigger, mwa_wl_energy)
    sites = map(sites) do s
        k = findall(t -> first(t_axis) + mwa_margin <= t <= last(t_axis) - mwa_margin, s.x)
        (x = s.x[k], x_tot = s.x_tot[k], max = s.max[k])
    end

    # output Table 
    return TypedTables.Table(
        blfc = blfc, timestamp = ts, eventID_fadc = evID, e_fc = efc,
        t_max = uconvert.(u"µs", estats.tmax), t_min = uconvert.(u"µs", estats.tmin), t_max_lar = uconvert.(u"µs", estats_trunc.tmax), t_min_lar = uconvert.(u"µs", estats_trunc.tmin),
        e_max = estats.max, e_min = estats.min, e_max_lar = estats_trunc.max, e_min_lar = estats_trunc.min,
        blmean = bl_stats.mean, blsigma = bl_stats.sigma, blslope = bl_stats.slope, bloffset = bl_stats.offset, 
        wfmean = sigstats.mean, wfsigma = sigstats.sigma, wfslope = sigstats.slope, wfoffset = sigstats.offset,
        # SG triggers
        threshold_sg = inters_thres, threshold_DC_sg = inters_thres_DC,
        trig_pos_sg = VectorOfVectors(inters.x), trig_max_sg = VectorOfVectors(inters.max),
        trig_pos_DC_sg = VectorOfVectors(inters_DC.x), trig_max_DC_sg = VectorOfVectors(inters_DC.max),
        # Trap triggers
        threshold_trap = inters_thres_trap, threshold_DC_trap = inters_thres_DC_trap,
        trig_pos_trap = VectorOfVectors(inters_trap.x), trig_pos_high_trap = VectorOfVectors(inters_trap.x_high),
        trig_pos_tot_trap = VectorOfVectors(inters_trap.x_tot), trig_max_trap = VectorOfVectors(inters_trap.max),
        trig_pos_DC_trap = VectorOfVectors(inters_DC_trap.x), trig_pos_high_DC_trap = VectorOfVectors(inters_DC_trap.x_high),
        trig_pos_tot_DC_trap = VectorOfVectors(inters_DC_trap.x_tot), trig_max_DC_trap = VectorOfVectors(inters_DC_trap.max),
        threshold_xt_tag = inters_thres_xt_tag,
        trig_pos_xt_tag = VectorOfVectors(inters_xt_tag.x), trig_pos_high_xt_tag = VectorOfVectors(inters_xt_tag.x_high),
        trig_pos_tot_xt_tag = VectorOfVectors(inters_xt_tag.x_tot), trig_max_xt_tag = VectorOfVectors(inters_xt_tag.max),
        # MWA triggers
        threshold_mwa = mwa_thresholds, trig_pos_mwa = VectorOfVectors([s.x for s in sites]),
        trig_pos_tot_mwa = VectorOfVectors([s.x_tot for s in sites]), trig_max_mwa = VectorOfVectors([s.max for s in sites])
    )
end
export dsp_sipm



"""
    dsp_sipm_compressed(data::Table, config::PropDict, pars::PropDict)

Variant of [`dsp_sipm`](@ref) for compressed data. The waveform is read from the column `waveform_bit_drop` instead of `waveform`; the arguments, the processing and the output columns, including the cross-talk tag columns, are the same as for `dsp_sipm`.

# Input columns
`waveform_bit_drop`, `baseline`, `timestamp`, `eventnumber`, `daqenergy`
"""
function dsp_sipm_compressed(data::Q, config::PropDict, pars::PropDict) where {Q <: Table}
    # common config
    sg_flt_degree         = config.sg_flt_degree
    t0_hpge_window        = config.t0_hpge_window

    # get SG config parameters
    sg_config = config.filters.sg
    min_tot_intersect     = sg_config.min_tot_intersect
    max_tot_intersect     = sg_config.max_tot_intersect
    min_threshold         = sg_config.min_threshold
    max_threshold         = sg_config.max_threshold
    n_σ_threshold         = sg_config.n_σ_threshold
    min_dc_threshold      = sg_config.min_dc_threshold
    max_dc_threshold      = sg_config.max_dc_threshold
    n_σ_dc_threshold      = sg_config.n_σ_dc_threshold

    # get trap config parameters
    trap_config = config.filters.trap
    trap_rt                = trap_config.rt
    trap_ft                = trap_config.ft
    trap_pz_tau            = trap_config.pz_tau
    trap_min_tot_intersect = trap_config.min_tot_intersect
    trap_max_tot_intersect = trap_config.max_tot_intersect
    trap_min_threshold     = trap_config.min_threshold
    trap_max_threshold     = trap_config.max_threshold
    trap_n_σ_threshold     = trap_config.n_σ_threshold
    trap_min_dc_threshold  = trap_config.min_dc_threshold
    trap_max_dc_threshold  = trap_config.max_dc_threshold
    trap_n_σ_dc_threshold  = trap_config.n_σ_dc_threshold

    xt_tag_config = config.filters.xt_tag
    xt_tag_window_length     = xt_tag_config.window_length
    xt_tag_n_σ_threshold     = xt_tag_config.n_σ_threshold
    xt_tag_min_tot_intersect = xt_tag_config.min_tot_intersect
    xt_tag_max_tot_intersect = xt_tag_config.max_tot_intersect

    # get MWA config parameters
    mwa_config = config.filters.mwa
    mwa_wl_trigger         = mwa_config.wl_trigger
    mwa_passes             = mwa_config.passes
    mwa_n_σ_threshold      = mwa_config.n_σ_threshold
    mwa_min_threshold      = mwa_config.min_threshold
    mwa_max_threshold      = mwa_config.max_threshold
    mwa_solve_overlap      = mwa_config.solve_overlap
    mwa_merge_distance     = mwa_config.merge_distance
    iseven(mwa_passes) || throw(ArgumentError("The MWA chain needs an even number of passes to stay centred, got $mwa_passes"))

    # unpack optimization parameters and, if the template step has run, the inverse kernel
    sg_window_length = pars.sg.wl
    mwa_wl_energy    = pars.mwa.wl
    has_kernel       = haskey(pars, :kernel)

    # get waveform data
    data = decode_data(data)
    wvfs = data.waveform_bit_drop
    blfc = data.baseline
    ts   = data.timestamp
    evID = data.eventnumber
    efc  = data.daqenergy

    # shift waveform by 0 to get Float64 conversion
    wvfs = shift_waveform.(wvfs, 0.0)

    # get wvf maximum and minimum with timepoints
    estats = extremestats.(wvfs)

    # get wvf maximum and minimum with timepoints for truncated waveform
    uflt_trunc = TruncateFilter(first(t0_hpge_window)..last(t0_hpge_window))
    estats_trunc = extremestats.(uflt_trunc.(wvfs))

    # === SG pipeline ===
    # savitzky golay filter: takes derivative of waveform plus smoothing
    sgflt_savitz = SavitzkyGolayFilter(sg_window_length, sg_flt_degree, 1)
    wvfs_sg = sgflt_savitz.(wvfs)

    # trigger finding on SG-filtered waveforms
    intflt_sg = IntersectMaximum(min_tot_intersect, max_tot_intersect)
    inters_thres = thresholdstats_mad.(wvfs_sg, min_threshold, max_threshold)
    inters = intflt_sg.(wvfs_sg, n_σ_threshold .* inters_thres)

    # integrate derivative
    integrator_filter = IntegratorFilter(gain=1)
    wvfs_int = integrator_filter.(wvfs_sg)

    # get blstats on integrated waveforms
    time_min = minimum(wvfs_int[1].time)
    time_max = last(wvfs_int[1].time)
    # baseline window ends at the first trigger (end of waveform if none), but spans at least 3 samples
    t_bl_min = time_min + 3*step(wvfs_int[1].time)
    t_bl_stop = map(x -> oftype(time_max, max(minimum(x; init = time_max), t_bl_min)), inters.x)
    bl_stats = signalstats.(wvfs_int, time_min, t_bl_stop)
    sigstats = signalstats.(wvfs_int, time_min, time_max)
    
    # flip around x-axis for discharge detection, relative to the baseline (median) of each integrated waveform
    wvfs_int_flip = multiply_waveform.(shift_waveform.(wvfs_int, .-median.(wvfs_int.signal)), -1.0)
    inters_thres_DC = thresholdstats_mad.(wvfs_int_flip, min_dc_threshold, max_dc_threshold)
    inters_DC = intflt_sg.(wvfs_int_flip, n_σ_dc_threshold .* inters_thres_DC)

    # === Trap pipeline ===
    # PZ correction on integrated waveforms
    pz_filter = InvCRFilter(trap_pz_tau)
    wvfs_pz = pz_filter.(wvfs_int)

    # trapezoidal charge filter
    trap_filter = TrapezoidalChargeFilter(trap_rt, trap_ft)
    wvfs_trap = trap_filter.(wvfs_pz)

    # trigger finding on trap-filtered waveforms
    intflt_trap = IntersectMaximum(trap_min_tot_intersect, trap_max_tot_intersect)
    inters_thres_trap = thresholdstats_mad.(wvfs_trap, trap_min_threshold, trap_max_threshold)
    inters_trap = intflt_trap.(wvfs_trap, trap_n_σ_threshold .* inters_thres_trap)

    # discharge detection on flipped integrated waveforms 
    inters_thres_DC_trap = thresholdstats_mad.(wvfs_int_flip, trap_min_dc_threshold, trap_max_dc_threshold)
    inters_DC_trap = intflt_trap.(wvfs_int_flip, trap_n_σ_dc_threshold .* inters_thres_DC_trap)

    wvfs_xt_tag = integrator_filter.(SavitzkyGolayFilter(xt_tag_window_length, sg_flt_degree, 1).(wvfs))
    wvfs_xt_tag = shift_waveform.(wvfs_xt_tag, .-median.(wvfs_xt_tag.signal))
    intflt_xt_tag = IntersectMaximum(xt_tag_min_tot_intersect, xt_tag_max_tot_intersect)
    inters_thres_xt_tag = thresholdstats_mad.(wvfs_xt_tag)
    inters_xt_tag = intflt_xt_tag.(wvfs_xt_tag, xt_tag_n_σ_threshold .* inters_thres_xt_tag)

    # === MWA pipeline ===
    # photon trace: the baseline-subtracted waveform through the inverse kernel of the template step; a channel without
    # template pars takes the derivative (kernel [1, -1]) instead
    t_axis = first(wvfs).time
    dt = step(t_axis)
    !has_kernel || pars.kernel.dt ≈ dt || throw(ArgumentError("The kernel taps are sampled at $(pars.kernel.dt), the waveforms at $dt"))
    photon_filter = CenteredConvolutionFilter(has_kernel ? Float64.(pars.kernel.taps) : [0.0, 1.0, -1.0])
    wvfs_bl = shift_waveform.(wvfs, [-median(wf.signal) for wf in wvfs])
    wvfs_photons = photon_filter.(wvfs_bl)

    # short (trigger) and long (energy) shaping, each scaled so that a unit stick (one photon after the kernel) has peak 1
    mwa_trigger_filter = MovingWindowMultiFilter(mwa_wl_trigger; passes = mwa_passes)
    mwa_energy_filter  = MovingWindowMultiFilter(mwa_wl_energy; passes = mwa_passes)
    stick = RDWaveform(t_axis, [i == length(t_axis) ÷ 2 ? 1.0 : 0.0 for i in eachindex(t_axis)])
    response_trigger = mwa_trigger_filter(stick).signal
    response_energy  = mwa_energy_filter(stick).signal
    wvfs_trigger = multiply_waveform.(mwa_trigger_filter.(wvfs_photons), 1 / maximum(response_trigger))
    wvfs_energy  = multiply_waveform.(mwa_energy_filter.(wvfs_photons), 1 / maximum(response_energy))
    wvfs_trigger = shift_waveform.(wvfs_trigger, [-median(wf.signal) for wf in wvfs_trigger])
    wvfs_energy  = shift_waveform.(wvfs_energy, [-median(wf.signal) for wf in wvfs_energy])

    # sites: the extrema beyond n σ of the trigger trace (matrix method) or of the energy trace
    wvfs_sites = mwa_solve_overlap ? wvfs_trigger : wvfs_energy
    mwa_thresholds = thresholdstats_mad.(wvfs_sites, mwa_min_threshold, mwa_max_threshold)
    sites = extrema_over_threshold.(wvfs_sites, mwa_n_σ_threshold .* mwa_thresholds)

    # matrix method: the amplitudes of all sites solved on the long shaping (otherwise the extrema of the energy trace)
    if mwa_solve_overlap
        solver = SiteSolver(response_energy ./ maximum(response_energy), argmax(response_energy), mwa_merge_distance * uconvert(NoUnits, mwa_wl_energy / dt))
        sites = map(zip(wvfs_energy, sites)) do (wf, s)
            sol = solver(wf.signal, uconvert.(NoUnits, (s.x .- first(t_axis)) ./ dt) .+ 1, uconvert.(NoUnits, s.x_tot ./ dt), Int.(sign.(s.max)))
            (x = convert.(eltype(t_axis), first(t_axis) .+ (sol.pos .- 1) .* dt), x_tot = convert.(eltype(t_axis), sol.tot .* dt), max = sol.A)
        end
    end

    # without the shaping transients at both ends of the trace
    mwa_margin = 2 * mwa_passes * max(mwa_wl_trigger, mwa_wl_energy)
    sites = map(sites) do s
        k = findall(t -> first(t_axis) + mwa_margin <= t <= last(t_axis) - mwa_margin, s.x)
        (x = s.x[k], x_tot = s.x_tot[k], max = s.max[k])
    end

    # output Table 
    return TypedTables.Table(
        blfc = blfc, timestamp = ts, eventID_fadc = evID, e_fc = efc,
        t_max = uconvert.(u"µs", estats.tmax), t_min = uconvert.(u"µs", estats.tmin), t_max_lar = uconvert.(u"µs", estats_trunc.tmax), t_min_lar = uconvert.(u"µs", estats_trunc.tmin),
        e_max = estats.max, e_min = estats.min, e_max_lar = estats_trunc.max, e_min_lar = estats_trunc.min,
        blmean = bl_stats.mean, blsigma = bl_stats.sigma, blslope = bl_stats.slope, bloffset = bl_stats.offset, 
        wfmean = sigstats.mean, wfsigma = sigstats.sigma, wfslope = sigstats.slope, wfoffset = sigstats.offset,
        # SG triggers
        threshold_sg = inters_thres, threshold_DC_sg = inters_thres_DC,
        trig_pos_sg = VectorOfVectors(inters.x), trig_max_sg = VectorOfVectors(inters.max),
        trig_pos_DC_sg = VectorOfVectors(inters_DC.x), trig_max_DC_sg = VectorOfVectors(inters_DC.max),
        # Trap triggers
        threshold_trap = inters_thres_trap, threshold_DC_trap = inters_thres_DC_trap,
        trig_pos_trap = VectorOfVectors(inters_trap.x), trig_pos_high_trap = VectorOfVectors(inters_trap.x_high),
        trig_pos_tot_trap = VectorOfVectors(inters_trap.x_tot), trig_max_trap = VectorOfVectors(inters_trap.max),
        trig_pos_DC_trap = VectorOfVectors(inters_DC_trap.x), trig_pos_high_DC_trap = VectorOfVectors(inters_DC_trap.x_high),
        trig_pos_tot_DC_trap = VectorOfVectors(inters_DC_trap.x_tot), trig_max_DC_trap = VectorOfVectors(inters_DC_trap.max),
        threshold_xt_tag = inters_thres_xt_tag,
        trig_pos_xt_tag = VectorOfVectors(inters_xt_tag.x), trig_pos_high_xt_tag = VectorOfVectors(inters_xt_tag.x_high),
        trig_pos_tot_xt_tag = VectorOfVectors(inters_xt_tag.x_tot), trig_max_xt_tag = VectorOfVectors(inters_xt_tag.max),
        # MWA triggers
        threshold_mwa = mwa_thresholds, trig_pos_mwa = VectorOfVectors([s.x for s in sites]),
        trig_pos_tot_mwa = VectorOfVectors([s.x_tot for s in sites]), trig_max_mwa = VectorOfVectors([s.max for s in sites])
    )
end
export dsp_sipm_compressed
