# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    get_t0(wvfs_pz::ArrayOfRDWaveforms, t0_threshold::Real; flt_pars::Vector{<:Unitful.Time}=[40u"ns", 100u"ns", 2000u"ns"], mintot::Unitful.Time=1500u"ns")

Start time t0 of the charge rise of each pole-zero-corrected waveform in `wvfs_pz`.

Each waveform is filtered with an asymmetric `TrapezoidalChargeFilter(flt_pars...)`. t0 is the interpolated time at which the filtered waveform first rises above `t0_threshold` and stays above it for at least `mintot`.

# Keywords
- `flt_pars`: pre-rise averaging time, gap time and post-rise averaging time of the trapezoidal filter
- `mintot`: minimum time over threshold

# Returns
- vector of t0 with unit µs, one entry per waveform. The entry is zero if the threshold is not crossed.
"""
function get_t0(wvfs_pz::ArrayOfRDWaveforms, t0_threshold::T; flt_pars::Vector{<:Unitful.Time{<:Real}}=[40u"ns", 100u"ns", 2000u"ns"], mintot::Unitful.Time{<:Real}=1500u"ns") where T<:Real
    # t0 determination
    # filter with fast asymetric trapezoidal filter and truncate waveform
    uflt_asy_t0 = TrapezoidalChargeFilter(flt_pars...)

    # eventuell zwei schritte!!!
    wvfs_flt_asy_t0 = uflt_asy_t0.(wvfs_pz)

    # get intersect at t0 threshold (fixed as in MJD analysis)
    flt_intersec_t0 = Intersect(mintot = mintot)

    # get t0 for every waveform as pick-off at fixed threshold
    t0 = uconvert.(u"µs", flt_intersec_t0.(wvfs_flt_asy_t0, t0_threshold).x)

    # replace NaN to avoid problems with further processing
    replace!(t0, NaN*unit(t0[1]) => zero(t0[1]))
end
export get_t0

"""
    get_threshold(wvfs::ArrayOfRDWaveforms, threshold::Union{Array{T}, T}; mintot::Unitful.Time=1000.0u"ns") where T<:Real

Time at which each waveform in `wvfs` first rises above `threshold`.

`threshold` is either one value for all waveforms or an array with one value per waveform. The time is interpolated linearly, and the waveform has to stay above the threshold for at least `mintot`.

# Returns
- vector of times with unit µs, one entry per waveform. The entry is zero if the threshold is not crossed.
"""
function get_threshold(wvfs::ArrayOfRDWaveforms, threshold::Union{Array{T}, T}; mintot::Unitful.Time=1000.0u"ns") where T<:Real
    # create intersect filter
    flt_intersect = Intersect(mintot = mintot)

    # get t0_threshold for every waveform as pick-off at fixed threshold
    t_thres = uconvert.(u"µs", flt_intersect.(wvfs, threshold).x)

    # replace NaN to avoid problems with further processing
    replace!(t_thres, NaN*unit(t_thres[1]) => zero(t_thres[1]))
end
export get_threshold

"""
    get_a(wvfs::ArrayOfRDWaveforms, current_window::ClosedInterval{<:Unitful.Time{<:Real}})

Maximum of the current waveforms `wvfs` in `current_window`, quadratically interpolated around the maximum sample
(`get_wvf_maximum`). Returns the maxima `max` and their timepoints `t`.
"""
get_a(wvfs::ArrayOfRDWaveforms, current_window::ClosedInterval{<:Unitful.Time{<:Real}}) = get_wvf_maximum.(wvfs, leftendpoint(current_window), rightendpoint(current_window))
export get_a


"""
    get_qdrift(wvfs::ArrayOfRDWaveforms, t_start::Array{<:Unitful.Time}, Δt::StepRangeLen{<:Unitful.Time}; pol_power::Int=3, sign_est_length::Unitful.Time=100u"ns")

Q-drift parameter of each waveform in `wvfs`: the difference of two charges collected after `t_start`.

The waveforms are integrated with `IntegratorFilter(1)`. Values of the integrated waveform at `t_start`, `t_start + first(Δt)` and `t_start + last(Δt)` are determined with a `SignalEstimator(PolynomialDNI(pol_power, sign_est_length))`. The result is the charge collected in `[t_start + first(Δt), t_start + last(Δt)]` minus the charge collected in `[t_start, t_start + first(Δt)]`.

# Arguments
- `t_start`: start time for each waveform, typically t0
- `Δt`: only the first and the last element are used, as the two integration times relative to `t_start`

# Keywords
- `pol_power`: polynomial degree of the signal estimator
- `sign_est_length`: length of the signal estimator window
"""
function get_qdrift(wvfs::ArrayOfRDWaveforms, t_start::Array{<:Unitful.Time{<:Real}}, Δt::StepRangeLen{<:Unitful.Time{<:Real}}; pol_power::Int=3, sign_est_length::Unitful.Time=100u"ns")
    # Integrate waveforms with gain = 1
    wvfs_flt_int = IntegratorFilter(1).(wvfs)
    
    # define signal estimator filter
    signal_estimator = SignalEstimator(PolynomialDNI(pol_power, sign_est_length))
    
    # calculate area between t_start and t_start + Δt
    area1  = signal_estimator.(wvfs_flt_int, t_start .+ first(Δt)) .- signal_estimator.(wvfs_flt_int, t_start)
    area2  = signal_estimator.(wvfs_flt_int, t_start .+ last(Δt)) .- signal_estimator.(wvfs_flt_int, t_start .+ first(Δt))
    
    # return qdrift as differenc of areas
    area2 .- area1
end
export get_qdrift

"""
    get_lq(wvfs_int::ArrayOfRDWaveforms, t_start::Array{<:Unitful.Time}, early::Unitful.Time, gap::Unitful.Time, late::Unitful.Time; pol_power::Int=3, sign_est_length::Unitful.Time=100u"ns")

LQ parameter: charge in the late window `[t_start + early + gap, t_start + early + gap + late]` minus the charge in the early window `[t_start, t_start + early]`.

`wvfs_int` are the integrated waveforms (`IntegratorFilter(1).(wvfs)`), so that the integration is done once for all LQ variants. The charges are determined with a `SignalEstimator(PolynomialDNI(pol_power, sign_est_length))`.

`get_qdrift(wvfs, t, Δt)` equals `get_lq(wvfs_int, t, first(Δt), 0u"µs", last(Δt) - first(Δt))`.
"""
function get_lq(wvfs_int::ArrayOfRDWaveforms, t_start::Array{<:Unitful.Time{<:Real}}, early::Unitful.Time, gap::Unitful.Time, late::Unitful.Time; pol_power::Int=3, sign_est_length::Unitful.Time=100u"ns")
    signal_estimator = SignalEstimator(PolynomialDNI(pol_power, sign_est_length))
    t_late = t_start .+ (early + gap)
    (signal_estimator.(wvfs_int, t_late .+ late) .- signal_estimator.(wvfs_int, t_late)) .- (signal_estimator.(wvfs_int, t_start .+ early) .- signal_estimator.(wvfs_int, t_start))
end
export get_lq

"""
    get_triggers_reversed(wvfs::ArrayOfRDWaveforms, sigma_threshold::Real, bl_window::ClosedInterval{<:Unitful.Time}; mintot::Unitful.Time=100.0u"ns")

Trigger positions and multiplicities found on the time-reversed waveforms.

The threshold of each waveform is `sigma_threshold` times its standard deviation in the baseline window, computed on the time axis of the first waveform from `first(wvfs[1].time) + leftendpoint(bl_window)` to `rightendpoint(bl_window)`. A threshold of zero is replaced by one. The waveforms are reversed in time and searched with `Intersect(mintot)`.

# Returns
- named tuple `(intersect, n)`: `intersect` holds the time of the first crossing in the original time axis, `n` the number of crossings of at least `mintot`.
"""
function get_triggers_reversed(wvfs::ArrayOfRDWaveforms, sigma_threshold::Real, bl_window::ClosedInterval{<:Unitful.Time{<:Real}}; mintot::Unitful.Time=100.0u"ns")
    flt_intersec_inTrace = Intersect(mintot = mintot)
    thres = signalstats.(wvfs, leftendpoint(bl_window) + first(wvfs[1].time), rightendpoint(bl_window)).sigma .* sigma_threshold
    # replace with high threshold if all values are saturated
    replace!(thres, zero(thres[1]) => one(thres[1]))

    inTrace_pileUp = flt_intersec_inTrace.(reverse_waveform.(wvfs), thres)

    (intersect = last(wvfs[1].time) .- inTrace_pileUp.x, n = inTrace_pileUp.multiplicity)
end
export get_triggers_reversed

"""
    get_triggers(wvfs::ArrayOfRDWaveforms, sigma_threshold::Real, bl_window::ClosedInterval{<:Unitful.Time}; mintot::Unitful.Time=100.0u"ns")

Trigger positions and multiplicities found on the waveforms.

The threshold of each waveform is `sigma_threshold` times its standard deviation in the baseline window, computed on the time axis of the first waveform from `first(wvfs[1].time) + leftendpoint(bl_window)` to `rightendpoint(bl_window)`. A threshold of zero is replaced by one. The waveforms are searched with `Intersect(mintot)`.

# Returns
- named tuple `(intersect, n)`: `intersect` holds the time of the first crossing (zero if there is none), `n` the number of crossings of at least `mintot`.
"""
function get_triggers(wvfs::ArrayOfRDWaveforms, sigma_threshold::Real, bl_window::ClosedInterval{<:Unitful.Time{<:Real}}; mintot::Unitful.Time=100.0u"ns")
    flt_intersec_inTrace = Intersect(mintot = mintot)
    thres = signalstats.(wvfs, leftendpoint(bl_window) + first(wvfs[1].time), rightendpoint(bl_window)).sigma .* sigma_threshold
    # replace with high threshold if all values are saturated
    replace!(thres, zero(thres[1]) => one(thres[1]))

    inTrace_pileUp = flt_intersec_inTrace.(wvfs, thres)

    (intersect = inTrace_pileUp.x, n = inTrace_pileUp.multiplicity)
end
export get_triggers
