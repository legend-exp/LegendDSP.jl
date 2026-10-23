# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    dsp_puls(data::Table, config::DSPConfig)

DSP chain for pulser waveforms. Returns a table with one row per waveform.

The baseline is the fit in `config.bl_window`. It is subtracted before the maximum, the time at 50% of the maximum and the energy of a trapezoidal filter are determined.

# Input columns
`waveform`, `baseline`, `timestamp`, `eventnumber`, `daqenergy`

# Output columns
- `blmean`, `blsigma`, `blslope`, `bloffset`: baseline statistics in `config.bl_window`
- `t50`: time at which the waveform reaches 50% of its maximum, in µs
- `e_max`: maximum of the baseline-subtracted waveform
- `e_10410`: maximum after a `TrapezoidalChargeFilter(10u"µs", 4u"µs")` (rise time 10 µs, flat-top time 4 µs)
- `blfc`: baseline from the FADC
- `timestamp`: timestamp from the FADC
- `eventID_fadc`: event number from the FADC
- `e_fc`: energy from the FADC
"""
function dsp_puls(data::Q, config::DSPConfig) where {Q <: Table}
    # get config parameters
    bl_window = config.bl_window

    # get waveform data
    data = decode_data(data)
    wvfs = data.waveform
    blfc = data.baseline
    ts   = data.timestamp
    evID = data.eventnumber
    efc  = data.daqenergy

    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs, leftendpoint(bl_window), rightendpoint(bl_window))

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)

    # get wvf maximum
    wvf_max = maximum.(wvfs.signal)

    # t50 determination
    t50 = get_threshold(wvfs, 0.5 .* wvf_max)

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





"""
    dsp_puls_compressed(data::Table, config::DSPConfig)

Variant of [`dsp_puls`](@ref) for compressed data. The waveform is read from the column `waveform_presummed` instead of `waveform`; the other input and all output columns are the same.
"""
function dsp_puls_compressed(data::Q, config::DSPConfig) where {Q <: Table}
    # get config parameters
    bl_window = config.bl_window

    # get waveform data
    data = decode_data(data)
    wvfs = data.waveform_presummed
    blfc = data.baseline
    ts   = data.timestamp
    evID = data.eventnumber
    efc  = data.daqenergy

    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs, leftendpoint(bl_window), rightendpoint(bl_window))

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)

    # get wvf maximum
    wvf_max = maximum.(wvfs.signal)

    # t50 determination
    t50 = get_threshold(wvfs, 0.5 .* wvf_max)

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
export dsp_puls_compressed


"""
    dsp_aux(data::Table, config::DSPConfig)

DSP chain for auxiliary channels whose waveforms are spikes. Returns a table with one row per waveform.

The `config` argument is accepted for consistency with the other DSP chains but is not used.

# Input columns
`waveform`, `baseline`, `timestamp`, `eventnumber`, `daqenergy`

# Output columns
- `e_max`: maximum sample of the waveform minus the FADC baseline
- `t_max`: time of the maximum sample
- `blfc`: baseline from the FADC
- `timestamp`: timestamp from the FADC
- `eventID_fadc`: event number from the FADC
- `e_fc`: energy from the FADC
"""
function dsp_aux(data::Q, ::DSPConfig) where {Q <: Table}
    data = decode_data(data)
    extrema = extremestats.(data.waveform)

    e_max = extrema.max .- data.baseline

    TypedTables.Table(
        e_max = e_max,
        t_max = extrema.tmax,
        blfc = data.baseline,
        timestamp = data.timestamp,
        eventID_fadc = data.eventnumber,
        e_fc = data.daqenergy,
    )
end
export dsp_aux


"""
    dsp_aux_compressed(data::Table, config::DSPConfig)

Variant of [`dsp_aux`](@ref) for compressed data. The waveform is read from the column `waveform_presummed`, and `e_max` is the maximum sample minus the FADC baseline multiplied by `presum_rate`. The input also needs the column `presum_rate`; the output columns are the same.
"""
function dsp_aux_compressed(data::Q, ::DSPConfig) where {Q <: Table}
    data = decode_data(data)
    wvfs = data.waveform_presummed
    extrema = extremestats.(wvfs)

    e_max = extrema.max .- data.baseline .* data.presum_rate

    TypedTables.Table(
        e_max = e_max,
        t_max = extrema.tmax,
        blfc = data.baseline,
        timestamp = data.timestamp,
        eventID_fadc = data.eventnumber,
        e_fc = data.daqenergy,
    )
end
export dsp_aux_compressed
