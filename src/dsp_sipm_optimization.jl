# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    dsp_sipm_thresholds_compressed(wvfs::ArrayOfRDWaveforms, config::PropDict)

This function calculates the baseline of the waveforms and the baseline of the waveforms with the sign flipped.
The function is used to calculate the thresholds for the SiPMs.

# Arguments
- `wvfs::ArrayOfRDWaveforms`: Array of RDWaveforms
- `config::PropDict`: Configuration parameters

# Returns
- `Table`: Table with the baseline of the waveforms and the baseline of the waveforms with the sign flipped
"""
function dsp_sipm_thresholds_compressed(wvfs::ArrayOfRDWaveforms, config::PropDict) where {Q <: Table}
    # get dsp meta parameters
    sg_window_length     = config.sg_window_length
    sg_flt_degree        = config.sg_flt_degree

    # get waveform data 
    wvfs = decode_data(wvfs)

    # shift waveform by 0 to get Float64 conversation --> ToDO: check if this is necessary
    wvfs = shift_waveform.(wvfs, 0.0)

    # savitzky golay filter: takes derivative of waveform plus smoothing
    sgflt_savitz = SavitzkyGolayFilter(sg_window_length, sg_flt_degree, 1)
    wvfs_sgflt_savitz = sgflt_savitz.(wvfs)

    # project waveforms on the y-axis
    bsl = vec(flatview(wvfs_sgflt_savitz.signal));

    # remove discharges
    # integrate derivative
    integrator_filter = IntegratorFilter(gain=1)
    wvfs_der_int = integrator_filter.(wvfs_sgflt_savitz)

    # flip around x-axis the filtered waveforms
    flipped_wf = multiply_waveform.(wvfs_der_int, -1.0)

    # find maxima in these flipped waveforms
    bsl_flipped = vec(flatview(flipped_wf.signal))

    # output Table 
    return TypedTables.Table(
        bsl = bsl,
        bsl_flipped = bsl_flipped
    )
end
export dsp_sipm_thresholds_compressed
