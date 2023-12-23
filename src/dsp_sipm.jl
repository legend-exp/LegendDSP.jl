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
    MINTOT_INTERSECT     = config.min_tot*u"ns"
    MAXTOT_INTERSECT     = config.max_tot*u"ns"
    SAVITZ_WINDOW_LENGTH = config.sg_window_length*u"ns"

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

    # savitzky golay filter: takes derivative of waveform plus smoothing
    sgflt_savitz = SavitzkyGolayFilter(SAVITZ_WINDOW_LENGTH, 2, 1)
    wvfs_sgflt_savitz = sgflt_savitz.(wvfs)

    # maximum finder
    intflt = IntersectMaximum(MINTOT_INTERSECT, MAXTOT_INTERSECT)
    inters = intflt.(wvfs_sgflt_savitz, threshold)

    # remove discharges
    # integrate and flip around x-axis the filtered waveforms
    integrator_filter = IntegratorFilter(gain=1)
    wvfs_der_int = integrator_filter.(wvfs_sgflt_savitz)
    flipped_wf = multiply_waveform.(wvfs_der_int, -1.0)

    inters_DC = intflt.(flipped_wf, threshold_DC)
    # filtered_inters = [inters[i] for i in 1:length(inters) if inters_DC[i].multiplicity == 0]
    # filtered_inters_struct = StructArray(filtered_inters) # has to be a struct array

    # output Table 
    return TypedTables.Table(
        blfc = blfc, timestamp = ts, eventID_fadc = evID, e_fc = efc,
        trig_pos =  VectorOfVectors(inters.x), trig_max =  VectorOfVectors(inters.max),
        trig_pos_DC =  VectorOfVectors(inters_DC.x), trig_max_DC =  VectorOfVectors(inters_DC.max)
    )
end
export dsp_sipm
