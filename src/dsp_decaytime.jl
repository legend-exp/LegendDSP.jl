# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    dsp_decay_times(wvfs::ArrayOfRDWaveforms, bl_window::ClosedInterval{<:Unitful.Time}, tail_window::ClosedInterval{<:Unitful.Time})
    dsp_decay_times(wvfs::ArrayOfRDWaveforms, config::DSPConfig)

Decay time of the exponential tail of each waveform in `wvfs`.

Each waveform is decoded, the mean in `bl_window` is subtracted, and [`tailstats`](@ref) fits a straight line to the logarithm of the samples in `tail_window`. The second method uses `config.bl_window` and `config.tail_window`.

# Returns
- vector of decay times τ with unit µs, one entry per waveform. The entry is zero if a sample in `tail_window` is not positive.
"""
function dsp_decay_times(wvfs::ArrayOfRDWaveforms, bl_window::ClosedInterval{<:Unitful.Time{<:T}}, tail_window::ClosedInterval{<:Unitful.Time{<:T}}) where T <: Real
    # get waveform data
    wvfs = decode_data(wvfs)

    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs, leftendpoint(bl_window), rightendpoint(bl_window))

    # substract baseline from waveforms
    wvfs_bl = shift_waveform.(wvfs, -bl_stats.mean)

    # extract decay times
    decay_times = tailstats.(wvfs_bl, leftendpoint(tail_window), rightendpoint(tail_window))
    
    # return converted to µs vals
    return uconvert.(u"µs", decay_times.τ)
end
export dsp_decay_times

dsp_decay_times(wvfs::ArrayOfRDWaveforms, config::DSPConfig) = dsp_decay_times(wvfs, config.bl_window, config.tail_window)