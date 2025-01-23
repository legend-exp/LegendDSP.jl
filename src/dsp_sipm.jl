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
function dsp_sipm(data::Q, config::PropDict, pars_threshold::PropDict) where {Q <: Table}
    # get dsp meta parameters
    min_tot_intersect     = config.min_tot_intersect
    max_tot_intersect     = config.max_tot_intersect
    sg_window_length     = config.sg_window_length
    sg_flt_degree        = config.sg_flt_degree
    t0_hpge_window       = config.t0_hpge_window

    # get config parameters
    threshold    = pars_threshold.sigma_thrs * config.n_sigma_threshold
    threshold_DC = pars_threshold.sigma_DC * config.n_sigma_dc_threshold

    # get waveform data 
    wvfs = data.waveform
    blfc = data.baseline
    ts   = data.timestamp
    evID = data.eventnumber
    efc  = data.daqenergy

    # shift waveform by 0 to get Float64 conversation --> ToDO: check if this is necessary
    wvfs = shift_waveform.(wvfs, 0.0)

    # get wvf maximum and minimum with timepoints
    flt_intersect = Intersect(mintot = 1u"ns")

    wvfs_max = maximum.(wvfs.signal)
    wvfs_min = minimum.(wvfs.signal)
    t_max = uconvert.(u"µs", flt_intersect.(wvfs, 0.99999 .* wvfs_max).x)
    t_min = uconvert.(u"µs", flt_intersect.(wvfs, wvfs_min).x)

    uflt_trunc = TruncateFilter(first(t0_hpge_window)..last(t0_hpge_window))
    wvfs_trunc = uflt_trunc.(wvfs)
    wvfs_trunc_max = maximum.(wvfs_trunc.signal)
    wvfs_trunc_min = minimum.(wvfs_trunc.signal)
    t_trunc_max = uconvert.(u"µs", flt_intersect.(wvfs_trunc, 0.99999 .* wvfs_trunc_max).x)
    t_trunc_min = uconvert.(u"µs", flt_intersect.(wvfs_trunc, wvfs_trunc_min).x)

    # savitzky golay filter: takes derivative of waveform plus smoothing
    sgflt_savitz = SavitzkyGolayFilter(sg_window_length, sg_flt_degree, 1)
    wvfs_sgflt_savitz = sgflt_savitz.(wvfs)

    # maximum finder
    intflt = IntersectMaximum(min_tot_intersect, max_tot_intersect)
    inters = intflt.(wvfs_sgflt_savitz, threshold)

    # remove discharges
    # integrate derivative
    integrator_filter = IntegratorFilter(gain=1)
    wvfs_der_int = integrator_filter.(wvfs_sgflt_savitz)

    # get blstats on derivative
    time_min = minimum(wvfs_der_int[1].time)
    Δt = 3*step(wvfs_der_int[1].time)
    bl_stats = signalstats.(wvfs_der_int, Ref(time_min), ifelse.(minimum.(inters.x; init=0u"s") .< time_min + Δt, time_min + Δt, minimum.(inters.x; init=0u"s")))
    sigstats = signalstats.(wvfs_der_int, time_min, last(wvfs_der_int[1].time))
    
    # flip around x-axis the filtered waveforms
    flipped_wf = multiply_waveform.(wvfs_der_int, -1.0)

    inters_DC = intflt.(flipped_wf, threshold_DC)

    # output Table 
    return TypedTables.Table(
        blfc = blfc, timestamp = ts, eventID_fadc = evID, e_fc = efc,
        t_max = t_max, t_min = t_min, t_max_lar = t_trunc_max, t_min_lar = t_trunc_min,
        e_max = wvfs_max, e_min = wvfs_min, e_max_lar = wvfs_trunc_max, e_min_lar = wvfs_trunc_min,
        blmean = bl_stats.mean, blsigma = bl_stats.sigma, blslope = bl_stats.slope, bloffset = bl_stats.offset, 
        wfmean = sigstats.mean, wfsigma = sigstats.sigma, wfslope = sigstats.slope, wfoffset = sigstats.offset,
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
function dsp_sipm_compressed(data::Q, config::PropDict, pars_threshold::PropDict) where {Q <: Table}
    # get dsp meta parameters
    min_tot_intersect     = config.min_tot_intersect
    max_tot_intersect     = config.max_tot_intersect
    sg_window_length = config.sg_window_length
    sg_flt_degree        = config.sg_flt_degree
    t0_hpge_window            = config.t0_hpge_window

    # get config parameters
    threshold    = pars_threshold.trig.σ * config.n_σ_threshold
    threshold_DC = pars_threshold.dc.σ * config.n_σ_dc_threshold

    # get waveform data 
    wvfs = decode_data(data.waveform_bit_drop)
    blfc = data.baseline
    ts   = data.timestamp
    evID = data.eventnumber
    efc  = data.daqenergy

    # shift waveform by 0 to get Float64 conversation --> ToDO: check if this is necessary
    wvfs = shift_waveform.(wvfs, 0.0)

    # get wvf maximum and minimum with timepoints
    flt_intersect = Intersect(mintot = 1u"ns")

    wvfs_max = maximum.(wvfs.signal)
    wvfs_min = minimum.(wvfs.signal)
    t_max = uconvert.(u"µs", flt_intersect.(wvfs, 0.99999 .* wvfs_max).x)
    t_min = uconvert.(u"µs", flt_intersect.(wvfs, wvfs_min).x)

    uflt_trunc = TruncateFilter(first(t0_hpge_window)..last(t0_hpge_window))
    wvfs_trunc = uflt_trunc.(wvfs)
    wvfs_trunc_max = maximum.(wvfs_trunc.signal)
    wvfs_trunc_min = minimum.(wvfs_trunc.signal)
    t_trunc_max = uconvert.(u"µs", flt_intersect.(wvfs_trunc, 0.99999 .* wvfs_trunc_max).x)
    t_trunc_min = uconvert.(u"µs", flt_intersect.(wvfs_trunc, wvfs_trunc_min).x)

    # savitzky golay filter: takes derivative of waveform plus smoothing
    sgflt_savitz = SavitzkyGolayFilter(sg_window_length, sg_flt_degree, 1)
    wvfs_sgflt_savitz = sgflt_savitz.(wvfs)

    # maximum finder
    intflt = IntersectMaximum(min_tot_intersect, max_tot_intersect)
    inters = intflt.(wvfs_sgflt_savitz, threshold)

    # remove discharges
    # integrate derivative
    integrator_filter = IntegratorFilter(gain=1)
    wvfs_der_int = integrator_filter.(wvfs_sgflt_savitz)

    # get blstats on derivative
    time_min = minimum(wvfs_der_int[1].time)
    Δt = 3*step(wvfs_der_int[1].time)
    bl_stats = signalstats.(wvfs_der_int, Ref(time_min), ifelse.(minimum.(inters.x; init=0u"s") .< time_min + Δt, time_min + Δt, minimum.(inters.x; init=0u"s")))
    sigstats = signalstats.(wvfs_der_int, time_min, last(wvfs_der_int[1].time))
    
    # flip around x-axis the filtered waveforms
    flipped_wf = multiply_waveform.(wvfs_der_int, -1.0)

    inters_DC = intflt.(flipped_wf, threshold_DC)

    # output Table 
    return TypedTables.Table(
        blfc = blfc, timestamp = ts, eventID_fadc = evID, e_fc = efc,
        t_max = t_max, t_min = t_min, t_max_lar = t_trunc_max, t_min_lar = t_trunc_min,
        e_max = wvfs_max, e_min = wvfs_min, e_max_lar = wvfs_trunc_max, e_min_lar = wvfs_trunc_min,
        blmean = bl_stats.mean, blsigma = bl_stats.sigma, blslope = bl_stats.slope, bloffset = bl_stats.offset, 
        wfmean = sigstats.mean, wfsigma = sigstats.sigma, wfslope = sigstats.slope, wfoffset = sigstats.offset,
        trig_pos =  VectorOfVectors(inters.x), trig_max =  VectorOfVectors(inters.max),
        trig_pos_DC =  VectorOfVectors(inters_DC.x), trig_max_DC =  VectorOfVectors(inters_DC.max)
    )
end
export dsp_sipm_compressed
