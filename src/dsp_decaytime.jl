# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    dsp_decay_times(wvfs::AbstractSamples, config::DSPConfig)

Get statistics on the logarhithmic of the tail of the `wvfs` in the interval `tail_window`.
# Returns
- `τ`: decay time in µs
"""
function dsp_decay_times(wvfs::ArrayOfRDWaveforms, config::DSPConfig)
    # get config parameters
    bl_window    = config.bl_window
    tail_window  = config.tail_window

    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs, first(bl_window), last(bl_window))

    # substract baseline from waveforms
    wvfs_bl = shift_waveform.(wvfs, -bl_stats.mean)

    # extract decay times
    decay_times = tailstats.(wvfs_bl, first(tail_window), last(tail_window))
    
    # return converted to µs vals
    return uconvert.(u"µs", decay_times.τ)
end
export dsp_decay_times