# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    dsp_sipm(data::Table, config::PropDict, pars_optimization::PropDict)

DSP chain for SiPM waveforms. Returns a table with one row per waveform.

# Input columns
`waveform`, `baseline`, `timestamp`, `eventnumber`, `daqenergy`

# Arguments
- `data`: table of waveforms; encoded columns are decoded with `decode_data`
- `config`: `PropDict` with the entries `sg_flt_degree`, `t0_hpge_window` (start and stop of the window in which the HPGe signal is expected, used for the `_lar` parameters) and `filters`. `filters.sg` holds the settings of the Savitzky-Golay pipeline and `filters.trap` those of the trapezoidal pipeline. Both contain `min_tot_intersect`, `max_tot_intersect`, `min_threshold`, `max_threshold`, `n_σ_threshold`, `min_dc_threshold`, `max_dc_threshold` and `n_σ_dc_threshold`; `filters.trap` also contains `rt`, `ft` and `pz_tau`.
- `pars_optimization`: `PropDict` with the optimized Savitzky-Golay window length `pars_optimization.sg.wl`

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

# Processing
- The waveform is converted to `Float64`.
- *Savitzky-Golay pipeline:* the waveform is filtered with a `SavitzkyGolayFilter(pars_optimization.sg.wl, config.sg_flt_degree, 1)`, which gives the smoothed derivative. The threshold is `n_σ_threshold` times [`thresholdstats_mad`](@ref) of the samples between `min_threshold` and `max_threshold`, and the triggers are found with [`IntersectMaximum`](@ref).
- The derivative is integrated, which gives a smoothed waveform, and the baseline statistics are determined on it.
- *Discharge detection:* the integrated waveform is multiplied by -1 and searched with the Savitzky-Golay `IntersectMaximum` using thresholds computed in the same way from the `dc` settings.
- *Trapezoidal pipeline:* the integrated waveform is pole-zero corrected with `pz_tau` and filtered with a `TrapezoidalChargeFilter(rt, ft)`. Thresholds and triggers are computed in the same way as for the Savitzky-Golay pipeline, with the settings of `filters.trap`.
- *Trapezoidal discharge detection:* the thresholds are computed from the flipped integrated waveform with the `dc` settings of `filters.trap`, and the triggers are found with the `IntersectMaximum` of the Savitzky-Golay settings.
"""
function dsp_sipm(data::Q, config::PropDict, pars_optimization::PropDict) where {Q <: Table}
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

    # unpack optimization parameters
    sg_window_length = pars_optimization.sg.wl

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
    
    # flip around x-axis for discharge detection
    wvfs_int_flip = multiply_waveform.(wvfs_int, -1.0)
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
    inters_DC_trap = intflt_sg.(wvfs_int_flip, trap_n_σ_dc_threshold .* inters_thres_DC_trap)

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
        trig_pos_tot_DC_trap = VectorOfVectors(inters_DC_trap.x_tot), trig_max_DC_trap = VectorOfVectors(inters_DC_trap.max)
    )
end
export dsp_sipm



"""
    dsp_sipm_compressed(data::Table, config::PropDict, pars_optimization::PropDict)

Variant of [`dsp_sipm`](@ref) for compressed data. The waveform is read from the column `waveform_bit_drop` instead of `waveform`; the arguments, the processing and the output columns are the same as for `dsp_sipm`.

# Input columns
`waveform_bit_drop`, `baseline`, `timestamp`, `eventnumber`, `daqenergy`
"""
function dsp_sipm_compressed(data::Q, config::PropDict, pars_optimization::PropDict) where {Q <: Table}
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

    # unpack optimization parameters
    sg_window_length = pars_optimization.sg.wl

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
    
    # flip around x-axis for discharge detection
    wvfs_int_flip = multiply_waveform.(wvfs_int, -1.0)
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
    inters_DC_trap = intflt_sg.(wvfs_int_flip, trap_n_σ_dc_threshold .* inters_thres_DC_trap)

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
        trig_pos_tot_DC_trap = VectorOfVectors(inters_DC_trap.x_tot), trig_max_DC_trap = VectorOfVectors(inters_DC_trap.max)
    )
end
export dsp_sipm_compressed
