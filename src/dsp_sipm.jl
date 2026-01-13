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
    # get dsp meta parameters
    min_tot_intersect     = config.min_tot_intersect
    max_tot_intersect     = config.max_tot_intersect
    sg_flt_degree         = config.sg_flt_degree
    t0_hpge_window        = config.t0_hpge_window
    min_threshold         = config.min_threshold
    max_threshold         = config.max_threshold
    n_σ_threshold         = config.n_σ_threshold
    min_dc_threshold      = config.min_dc_threshold
    max_dc_threshold      = config.max_dc_threshold
    n_σ_dc_threshold      = config.n_σ_dc_threshold

    # unpack optimization parameters - extract numerical value from PropDict
    sg_window_length = haskey(pars_optimization.sg.wl, :val) ? pars_optimization.sg.wl.val : pars_optimization.sg.wl

    # get waveform data 
    wvfs = data.waveform
    blfc = data.baseline
    ts   = data.timestamp
    evID = data.eventnumber
    efc  = data.daqenergy

    # shift waveform by 0 to get Float64 conversation --> ToDO: check if this is necessary
    wvfs = shift_waveform.(wvfs, 0.0)

    # get wvf maximum and minimum with timepoints
    estats = extremestats.(wvfs)

    # get wvf maximum and minimum with timepoints for truncated waveform
    uflt_trunc = TruncateFilter(first(t0_hpge_window)..last(t0_hpge_window))
    estats_trunc = extremestats.(uflt_trunc.(wvfs))

    # savitzky golay filter: takes derivative of waveform plus smoothing
    sgflt_savitz = SavitzkyGolayFilter(sg_window_length, sg_flt_degree, 1)
    wvfs = sgflt_savitz.(wvfs)

    # maximum finder
    intflt = IntersectMaximum(min_tot_intersect, max_tot_intersect)
    inters_thres = thresholdstats.(wvfs, min_threshold, max_threshold)
    inters = intflt.(wvfs, n_σ_threshold .* inters_thres)

    # remove discharges
    # integrate derivative
    integrator_filter = IntegratorFilter(gain=1)
    wvfs = integrator_filter.(wvfs)

    # get blstats on derivative
    time_min = minimum(wvfs[1].time)
    Δt = 3*step(wvfs[1].time)
    bl_stats = signalstats.(wvfs, Ref(time_min), ifelse.(minimum.(inters.x; init=0u"s") .< time_min + Δt, time_min + Δt, minimum.(inters.x; init=0u"s")))
    sigstats = signalstats.(wvfs, time_min, last(wvfs[1].time))
    
    # flip around x-axis the filtered waveforms
    wvfs = multiply_waveform.(wvfs, -1.0)
    inters_thres_DC = thresholdstats.(wvfs, min_dc_threshold, max_dc_threshold)
    inters_DC = intflt.(wvfs, n_σ_dc_threshold .* inters_thres_DC)

    # output Table 
    return TypedTables.Table(
        blfc = blfc, timestamp = ts, eventID_fadc = evID, e_fc = efc,
        t_max = uconvert.(u"µs", estats.tmax), t_min = uconvert.(u"µs", estats.tmin), t_max_lar = uconvert.(u"µs", estats_trunc.tmax), t_min_lar = uconvert.(u"µs", estats_trunc.tmin),
        e_max = wvfs_max, e_min = wvfs_min, e_max_lar = wvfs_trunc_max, e_min_lar = wvfs_trunc_min,
        blmean = bl_stats.mean, blsigma = bl_stats.sigma, blslope = bl_stats.slope, bloffset = bl_stats.offset, 
        wfmean = sigstats.mean, wfsigma = sigstats.sigma, wfslope = sigstats.slope, wfoffset = sigstats.offset,
        threshold = inters_thres, threshold_DC = inters_thres_DC,
        trig_pos =  VectorOfVectors(inters.x), trig_max =  VectorOfVectors(inters.max),
        trig_pos_DC =  VectorOfVectors(inters_DC.x), trig_max_DC =  VectorOfVectors(inters_DC.max)
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
    # get dsp meta parameters
    min_tot_intersect     = config.min_tot_intersect
    max_tot_intersect     = config.max_tot_intersect
    sg_flt_degree         = config.sg_flt_degree
    t0_hpge_window        = config.t0_hpge_window
    min_threshold         = config.min_threshold
    max_threshold         = config.max_threshold
    n_σ_threshold         = config.n_σ_threshold
    min_dc_threshold      = config.min_dc_threshold
    max_dc_threshold      = config.max_dc_threshold
    n_σ_dc_threshold      = config.n_σ_dc_threshold

    # unpack optimization parameters
    sg_window_length = pars_optimization.sg.wl

    # get waveform data 
    wvfs = decode_data(data.waveform_bit_drop)
    blfc = data.baseline
    ts   = data.timestamp
    evID = data.eventnumber
    efc  = data.daqenergy

    # shift waveform by 0 to get Float64 conversation --> ToDO: check if this is necessary
    wvfs = shift_waveform.(wvfs, 0.0)

    # get wvf maximum and minimum with timepoints
    estats = extremestats.(wvfs)

    # get wvf maximum and minimum with timepoints for truncated waveform
    uflt_trunc = TruncateFilter(first(t0_hpge_window)..last(t0_hpge_window))
    estats_trunc = extremestats.(uflt_trunc.(wvfs))

    # savitzky golay filter: takes derivative of waveform plus smoothing
    sgflt_savitz = SavitzkyGolayFilter(sg_window_length, sg_flt_degree, 1)
    wvfs = sgflt_savitz.(wvfs)

    # maximum finder
    intflt = IntersectMaximum(min_tot_intersect, max_tot_intersect)
    inters_thres = thresholdstats.(wvfs, min_threshold, max_threshold)
    inters = intflt.(wvfs, n_σ_threshold .* inters_thres)

    # remove discharges
    # integrate derivative
    integrator_filter = IntegratorFilter(gain=1)
    wvfs = integrator_filter.(wvfs)

    # get blstats on derivative
    time_min = minimum(wvfs[1].time)
    Δt = 3*step(wvfs[1].time)
    bl_stats = signalstats.(wvfs, Ref(time_min), ifelse.(minimum.(inters.x; init=0u"s") .< time_min + Δt, time_min + Δt, minimum.(inters.x; init=0u"s")))
    sigstats = signalstats.(wvfs, time_min, last(wvfs[1].time))
    
    # flip around x-axis the filtered waveforms
    wvfs = multiply_waveform.(wvfs, -1.0)
    inters_thres_DC = thresholdstats.(wvfs, min_dc_threshold, max_dc_threshold)
    inters_DC = intflt.(wvfs, n_σ_dc_threshold .* inters_thres_DC)

    # output Table 
    return TypedTables.Table(
        blfc = blfc, timestamp = ts, eventID_fadc = evID, e_fc = efc,
        t_max = uconvert.(u"µs", estats.tmax), t_min = uconvert.(u"µs", estats.tmin), t_max_lar = uconvert.(u"µs", estats_trunc.tmax), t_min_lar = uconvert.(u"µs", estats_trunc.tmin),
        e_max = estats.max, e_min = estats.min, e_max_lar = estats_trunc.max, e_min_lar = estats_trunc.min,
        blmean = bl_stats.mean, blsigma = bl_stats.sigma, blslope = bl_stats.slope, bloffset = bl_stats.offset, 
        wfmean = sigstats.mean, wfsigma = sigstats.sigma, wfslope = sigstats.slope, wfoffset = sigstats.offset,
        threshold = inters_thres, threshold_DC = inters_thres_DC,
        trig_pos =  VectorOfVectors(inters.x), trig_max =  VectorOfVectors(inters.max),
        trig_pos_DC =  VectorOfVectors(inters_DC.x), trig_max_DC =  VectorOfVectors(inters_DC.max)
    )
end
export dsp_sipm_compressed
