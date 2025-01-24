# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

function dsp_pmts(data::Q, config::PropDict) where {Q <: Table}
    # get dsp meta parameters
    baseline_window_start = config.default.baseline_window_start
    baseline_window_end   = config.default.baseline_window_end
    min_tot_intersect     = config.default.min_tot_intersect
    max_tot_intersect     = config.default.max_tot_intersect
    intersect_threshold   = config.default.intersect_threshold
    wsg_window_length     = config.default.sg_window_length
    wsg_flt_degree        = config.default.sg_flt_degree
    wsg_weight            = config.default.wsg_weight
    saturation_limit_high  = config.default.saturation_limit_high
    saturation_limit_low   = config.default.saturation_limit_low

    # get waveform data 
    wvfs = data.waveform
    ts   = data.timestamp
    ch  = data.channel
    blstats = signalstats.(wvfs, baseline_window_start, baseline_window_end) 

    # substract baseline_window_end
    wvfs = shift_waveform.(wvfs, -blstats.mean)

    # get wvf maximum and minimum with timepoints
    raw_pulse_params = extremestats.(wf_blsub)
    
    # get every peak above threshold
    intflt = IntersectMaximum(min_tot_intersect, max_tot_intersect) 
    trig = intflt.(wf_blsub, intersect_threshold) 

    # get saturated peaks
    sat = saturation.(wf_blsub, saturation_limit_low, saturation_limit_high) 

    # weighted savitzky golay filter
    w_sg = WeightedSavitzkyGolayFilter(wsg_window_length, wsg_flt_degree, wsg_weight) 
    w_sg_filter = w_sg.(wf_blsub)

    # get smoothed wvf maximum and minimum with timepoints
    pulse_params = extremestats.(w_sg_filter)

    # output Table 
    return TypedTables.Table(
        timestamp = ts, channel = ch, rawPulseHeight = raw_pulse_params.max, rawPulseLow = raw_pulse_params.min,
        raw_tp_hi = raw_pulse_params.tmax, raw_tp_low = raw_pulse_params.tmin,
        trig_max = trig.max, trig_t = trig.x, trig_mult = trig.multiplicity,
        sat_low = sat.low, sat_high = sat.high,
        pulseHeight = pulse_params.max, pulseLow = pulse_params.min,
        tp_hi = pulse_params.tmax, tp_low = pulse_params.tmin,
        bl_mean = bl_stats.mean, bl_sigma = bl_stats.sigma, bl_slope = bl_stats.slope
    )
end