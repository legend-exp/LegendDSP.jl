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
    uconvert.(u"µs", flt_intersec_t0.(wvfs_flt_asy_t0, t0_threshold).x)
end

"""
    get_t50(wvfs_pz::ArrayOfRDWaveforms, wvf_max::T) where T<:Real

Get t50 for each waveform in `wvfs_pz` by intersecting the waveform with a fixed threshold at `wvf_max * 0.5`.
"""
function get_t50(wvfs_pz::ArrayOfRDWaveforms, wvf_max::Array{T}) where T<:Real
    # create intersect filter
    flt_intersect = Intersect(mintot = 1000u"ns")

    # get t50 for every waveform as pick-off at fixed threshold
    uconvert.(u"µs", flt_intersect.(wvfs_pz, wvf_max .* 0.5).x)
end