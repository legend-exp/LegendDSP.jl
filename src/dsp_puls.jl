# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

""" 
    dsp_puls(data::Q, config::DSPConfig) where {Q <: Table}

DSP function for pulser processing.
# Input data
The input data is a table with the following columns:
- `waveform`: waveform data
- `baseline`: baseline data
- `timestamp`: timestamp data
- `eventnumber`: event number data
- `daqenergy`: energy data

# Output data
The output data is a table with the following columns:
- `blmean`: baseline mean
- `blsigma`: baseline sigma
- `blslope`: baseline slope
- `bloffset`: baseline offset
- `t50`: timepoint of 50% of waveform maximum
- `e_max`: maximum of waveform
- `e_10410`: energy of waveform with trapezoidal filter of 10µs rise time with 4µs flat-top
- `blfc`: baseline from FADC
- `timestamp`: timestamp
- `eventID_fadc`: event number from FADC
- `e_fc`: energy from FADC
"""
function dsp_puls(data::Q, config::DSPConfig) where {Q <: Table}
    # get config parameters
    bl_window = config.bl_window

    # get waveform data 
    wvfs = data.waveform
    blfc = data.baseline
    ts   = data.timestamp
    evID = data.eventnumber
    efc  = data.daqenergy

    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs, first(bl_window), last(bl_window))

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)

    # get wvf maximum
    wvf_max = maximum.(wvfs.signal)

    # t50 determination
    t50 = get_t50(wvfs, wvf_max)

    # extract energy and ENC noise param from maximum of filtered wvfs
    uflt_10410 = TrapezoidalChargeFilter(10u"µs", 4u"µs")

    wvfs_flt = uflt_10410.(wvfs)
    e_10410  = maximum.(wvfs_flt.signal)

    # output Table 
    return TypedTables.Table(blmean = bl_stats.mean, blsigma = bl_stats.sigma, blslope = bl_stats.slope, bloffset = bl_stats.offset, 
    t50 = t50,
    e_max = wvf_max, 
    e_10410 = e_10410, 
    blfc = blfc, timestamp = ts, eventID_fadc = evID, e_fc = efc,
    )
end
export dsp_puls