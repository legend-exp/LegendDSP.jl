# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

function dsp_pmts(data::Q, config::PropDict) where {Q <: Table}
    # get dsp meta parameters
    time_axis_step_length = config.time_axis_step_length
    baseline_window_start = config.baseline_window_start
    baseline_window_end   = config.baseline_window_end
    min_tot_intersect     = config.min_tot_intersect
    max_tot_intersect     = config.max_tot_intersect
    intersect_threshold   = config.intersect_threshold
    wsg_window_length     = config.wsg_window_length
    wsg_flt_degree        = config.wsg_flt_degree
    wsg_weight            = config.wsg_weight
    saturation_limit_high = config.saturation_limit_high
    saturation_limit_low  = config.saturation_limit_low

    # get waveform data 
    waveform = decode_data(data.waveform)
    ts  = data.timestamp
    ch  = data.channel
    evID = data.eventnumber
    efc  = data.daqenergy
    
    time_flt = TimeAxisFilter(time_axis_step_length)
    wvfs = time_flt.(waveform)
    
    # calculate baseline stats
    bl_stats = signalstats.(wvfs, baseline_window_start, baseline_window_end) 

    # substract baseline
    wf_blsub = shift_waveform.(wvfs, -bl_stats.mean)

    # get wvf maximum and minimum with timepoints
    raw_pulse_params = extremestats.(wf_blsub)
    
    # get every peak above threshold
    intflt = IntersectMaximum(min_tot_intersect, max_tot_intersect) 
    trig = intflt.(wf_blsub, intersect_threshold) 

    # get saturated peaks
    sat = saturation.(wvfs, saturation_limit_low, saturation_limit_high) 

    # weighted savitzky golay filter
    w_sg = if wsg_weight == 0
        SavitzkyGolayFilter(wsg_window_length, Int(wsg_flt_degree), 0)
    else
        WeightedSavitzkyGolayFilter(wsg_window_length, Int(wsg_flt_degree), Int(wsg_weight))
    end
    w_sg_filter = w_sg.(wf_blsub)

    # get smoothed wvf maximum and minimum with timepoints
    pulse_params = extremestats.(w_sg_filter)

    # output Table 
    return TypedTables.Table(
        timestamp = ts, eventID_fadc = evID, e_fc = efc, channel = ch,
        raw_pulse_height = raw_pulse_params.max, raw_pulse_low = raw_pulse_params.min,
        raw_t0_hi = raw_pulse_params.tmax, raw_t0_low = raw_pulse_params.tmin,
        trig_max = VectorOfVectors(trig.max), trig_t = VectorOfVectors(trig.x), trig_mult = trig.multiplicity,
        sat_low = sat.low, sat_high = sat.high,
        pulse_height = pulse_params.max, pulse_low = pulse_params.min,
        t0_hi = pulse_params.tmax, t0_low = pulse_params.tmin,
        bl_mean = bl_stats.mean, bl_sigma = bl_stats.sigma, bl_slope = bl_stats.slope
    )
end
export dsp_pmts