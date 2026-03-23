# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    dsp_sipm(data::Q, config::PropDict, pars_threshold::PropDict)

DSP routine for SiPM data. It needs the threshold parameters from the threshold scan for the given SiPM channel as well as the discharge threshold parameters.

# Input data
The input data is a table with the following columns:
- `waveform`: waveform data
- `baseline`: baseline data
- `timestamp`: timestamp data
- `eventnumber`: event number data
- `daqenergy`: energy data

# Output data
The output data is a table with the following columns:
- `blfc`: baseline from FADC
- `timestamp`: timestamp
- `eventID_fadc`: event number from FADC
- `e_fc`: energy from FADC
- `trig_pos`: trigger positions from DSP
- `trig_max`: trigger maxima from DSP
- `trig_pos_DC`: trigger positions of discharges
- `trig_max_DC`: trigger maxima of discharges
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
    sg_threshold_method   = get(sg_config, :threshold_method, "std")

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
    trap_threshold_method = get(trap_config, :threshold_method, "std")

    # unpack optimization parameters
    sg_window_length = pars_optimization.sg.wl

    # get waveform data 
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
    inters_thres = (sg_threshold_method == "mad" ? thresholdstats_mad : thresholdstats).(wvfs_sg, min_threshold, max_threshold)
    inters = intflt_sg.(wvfs_sg, n_σ_threshold .* inters_thres)

    # integrate derivative
    integrator_filter = IntegratorFilter(gain=1)
    wvfs_int = integrator_filter.(wvfs_sg)

    # get blstats on integrated waveforms
    time_min = minimum(wvfs_int[1].time)
    Δt = 3*step(wvfs_int[1].time)
    bl_stats = signalstats.(wvfs_int, Ref(time_min), ifelse.(minimum.(inters.x; init=0u"s") .< time_min + Δt, time_min + Δt, minimum.(inters.x; init=0u"s")))
    sigstats = signalstats.(wvfs_int, time_min, last(wvfs_int[1].time))
    
    # flip around x-axis for discharge detection
    wvfs_int_flip = multiply_waveform.(wvfs_int, -1.0)
    inters_thres_DC = (sg_threshold_method == "mad" ? thresholdstats_mad : thresholdstats).(wvfs_int_flip, min_dc_threshold, max_dc_threshold)
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
    inters_thres_trap = (trap_threshold_method == "mad" ? thresholdstats_mad : thresholdstats).(wvfs_trap, trap_min_threshold, trap_max_threshold)
    inters_trap = intflt_trap.(wvfs_trap, trap_n_σ_threshold .* inters_thres_trap)

    # discharge detection on flipped integrated waveforms 
    inters_thres_DC_trap = (trap_threshold_method == "mad" ? thresholdstats_mad : thresholdstats).(wvfs_int_flip, trap_min_dc_threshold, trap_max_dc_threshold)
    inters_DC_trap = intflt_sg.(wvfs_int_flip, trap_n_σ_dc_threshold .* inters_thres_DC_trap)

    # output Table 
    return TypedTables.Table(
        blfc = blfc, timestamp = ts, eventID_fadc = evID, e_fc = efc,
        t_max = uconvert.(u"µs", estats.tmax), t_min = uconvert.(u"µs", estats.tmin), t_max_lar = uconvert.(u"µs", estats_trunc.tmax), t_min_lar = uconvert.(u"µs", estats_trunc.tmin),
        e_max = estats.max, e_min = estats.min, e_max_lar = estats_trunc.max, e_min_lar = estats_trunc.min,
        blmean = bl_stats.mean, blsigma = bl_stats.sigma, blslope = bl_stats.slope, bloffset = bl_stats.offset, 
        wfmean = sigstats.mean, wfsigma = sigstats.sigma, wfslope = sigstats.slope, wfoffset = sigstats.offset,
        # SG triggers
        threshold = inters_thres, threshold_DC = inters_thres_DC,
        trig_pos = VectorOfVectors(inters.x), trig_max = VectorOfVectors(inters.max),
        trig_pos_DC = VectorOfVectors(inters_DC.x), trig_max_DC = VectorOfVectors(inters_DC.max),
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
    dsp_sipm_compressed(data::Q, config::PropDict, pars_threshold::PropDict)

DSP routine for SiPM data. It needs the threshold parameters from the threshold scan for the given SiPM channel as well as the discharge threshold parameters.

# Input data
The input data is a table with the following columns:
- `waveform`: waveform data
- `baseline`: baseline data
- `timestamp`: timestamp data
- `eventnumber`: event number data
- `daqenergy`: energy data

# Output data
The output data is a table with the following columns:
- `blfc`: baseline from FADC
- `timestamp`: timestamp
- `eventID_fadc`: event number from FADC
- `e_fc`: energy from FADC
- `trig_pos`: trigger positions from DSP
- `trig_max`: trigger maxima from DSP
- `trig_pos_DC`: trigger positions of discharges
- `trig_max_DC`: trigger maxima of discharges
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
    sg_threshold_method   = get(sg_config, :threshold_method, "std")

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
    trap_threshold_method = get(trap_config, :threshold_method, "std")

    # unpack optimization parameters
    sg_window_length = pars_optimization.sg.wl

    # get waveform data 
    wvfs = decode_data(data.waveform_bit_drop)
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
    inters_thres = (sg_threshold_method == "mad" ? thresholdstats_mad : thresholdstats).(wvfs_sg, min_threshold, max_threshold)
    inters = intflt_sg.(wvfs_sg, n_σ_threshold .* inters_thres)

    # integrate derivative
    integrator_filter = IntegratorFilter(gain=1)
    wvfs_int = integrator_filter.(wvfs_sg)

    # get blstats on integrated waveforms
    time_min = minimum(wvfs_int[1].time)
    Δt = 3*step(wvfs_int[1].time)
    bl_stats = signalstats.(wvfs_int, Ref(time_min), ifelse.(minimum.(inters.x; init=0u"s") .< time_min + Δt, time_min + Δt, minimum.(inters.x; init=0u"s")))
    sigstats = signalstats.(wvfs_int, time_min, last(wvfs_int[1].time))
    
    # flip around x-axis for discharge detection
    wvfs_int_flip = multiply_waveform.(wvfs_int, -1.0)
    inters_thres_DC = (sg_threshold_method == "mad" ? thresholdstats_mad : thresholdstats).(wvfs_int_flip, min_dc_threshold, max_dc_threshold)
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
    inters_thres_trap = (trap_threshold_method == "mad" ? thresholdstats_mad : thresholdstats).(wvfs_trap, trap_min_threshold, trap_max_threshold)
    inters_trap = intflt_trap.(wvfs_trap, trap_n_σ_threshold .* inters_thres_trap)

    # discharge detection on flipped integrated waveforms 
    inters_thres_DC_trap = (trap_threshold_method == "mad" ? thresholdstats_mad : thresholdstats).(wvfs_int_flip, trap_min_dc_threshold, trap_max_dc_threshold)
    inters_DC_trap = intflt_sg.(wvfs_int_flip, trap_n_σ_dc_threshold .* inters_thres_DC_trap)

    # output Table 
    return TypedTables.Table(
        blfc = blfc, timestamp = ts, eventID_fadc = evID, e_fc = efc,
        t_max = uconvert.(u"µs", estats.tmax), t_min = uconvert.(u"µs", estats.tmin), t_max_lar = uconvert.(u"µs", estats_trunc.tmax), t_min_lar = uconvert.(u"µs", estats_trunc.tmin),
        e_max = estats.max, e_min = estats.min, e_max_lar = estats_trunc.max, e_min_lar = estats_trunc.min,
        blmean = bl_stats.mean, blsigma = bl_stats.sigma, blslope = bl_stats.slope, bloffset = bl_stats.offset, 
        wfmean = sigstats.mean, wfsigma = sigstats.sigma, wfslope = sigstats.slope, wfoffset = sigstats.offset,
        # SG triggers
        threshold = inters_thres, threshold_DC = inters_thres_DC,
        trig_pos = VectorOfVectors(inters.x), trig_max = VectorOfVectors(inters.max),
        trig_pos_DC = VectorOfVectors(inters_DC.x), trig_max_DC = VectorOfVectors(inters_DC.max),
        # Trap triggers
        threshold_trap = inters_thres_trap, threshold_DC_trap = inters_thres_DC_trap,
        trig_pos_trap = VectorOfVectors(inters_trap.x), trig_pos_high_trap = VectorOfVectors(inters_trap.x_high),
        trig_pos_tot_trap = VectorOfVectors(inters_trap.x_tot), trig_max_trap = VectorOfVectors(inters_trap.max),
        trig_pos_DC_trap = VectorOfVectors(inters_DC_trap.x), trig_pos_high_DC_trap = VectorOfVectors(inters_DC_trap.x_high),
        trig_pos_tot_DC_trap = VectorOfVectors(inters_DC_trap.x_tot), trig_max_DC_trap = VectorOfVectors(inters_DC_trap.max)
    )
end
export dsp_sipm_compressed
