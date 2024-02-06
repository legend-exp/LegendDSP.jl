# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    get_t0(wvfs_pz::ArrayOfRDWaveforms, t0_threshold::T) where T<:Real

Get t0 for each waveform in `wvfs_pz` by using a fast asymetric trapezoidal filter and
a fixed threshold at `t0_threshold`. The filter is truncated to the range 0µs to 60µs where the Ge trigger is expected in FlashCam.
"""
function get_t0(wvfs_pz::ArrayOfRDWaveforms, t0_threshold::T) where T<:Real
    # t0 determination
    # filter with fast asymetric trapezoidal filter and truncate waveform
    uflt_asy_t0 = TrapezoidalChargeFilter(40u"ns", 100u"ns", 2000u"ns")
    uflt_trunc_t0 = TruncateFilter(0u"µs"..60u"µs")

    # eventuell zwei schritte!!!
    wvfs_flt_asy_t0 = uflt_asy_t0.(uflt_trunc_t0.(wvfs_pz))

    # get intersect at t0 threshold (fixed as in MJD analysis)
    flt_intersec_t0 = Intersect(mintot = 1500u"ns")

    # get t0 for every waveform as pick-off at fixed threshold
    t0 = uconvert.(u"µs", flt_intersec_t0.(wvfs_flt_asy_t0, t0_threshold).x)

    # replace NaN to avoid problems with further processing
    replace!(t0, NaN*unit(t0[1]) => zero(t0[1]))
end
export get_t0

"""
    get_threshold(wvfs::ArrayOfRDWaveforms, threshold::Array{T}) where T<:Real

Get threshold for each waveform in `wvfs` by intersecting the waveform with a `threshold` per waveform or globally.
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
    get_qdrift(wvfs::ArrayOfRDWaveforms, t_start::Array{Unitful.Time{T}}, Δt::UnitRange{Unitful.Time{T}}; pol_power::Int=3, sign_est_length::Unitful.Time=100u"ns")

Get the Q-drift parameter for each waveform in `wvfs` by integrating the waveform with gain = 1 and using a polynomial signal estimator of order `pol_power` and length `sign_est_length` to estimate the signal.
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
    get_intracePileUp(wvfs::ArrayOfRDWaveforms, sigma_threshold::Real, blmin::Unitful.Time{<:Real}, blmax::Unitful.Time{<:Real}; mintot::Unitful.Time=100.0u"ns")

Get position and multiplicity of in-trace pile-up as intersect of reversed derivative signal with threshold as multiple of std. The `wvfs` have to be a current signal.
"""
function get_intracePileUp(wvfs::ArrayOfRDWaveforms, sigma_threshold::Real, bl_window::ClosedInterval{Unitful.Time{<:Real}}; mintot::Unitful.Time=100.0u"ns")
    # get position and multiplicity of in-trace pile-up as intersect of reversed derivative signal with threshold as multiple of std
    flt_intersec_inTrace = Intersect(mintot = mintot)
    inTrace_pileUp       = flt_intersec_inTrace.(reverse_waveform.(wvfs), signalstats.(wvfs, first(bl_window) + first(wvfs[1].time), last(bl_window)).sigma .* sigma_threshold)
    # return intersect position measured from the non-reversed waveform and multiplicity
    (intersect = last(wvfs[1].time) .- inTrace_pileUp.x, n = inTrace_pileUp.multiplicity)
end
export get_intracePileUp