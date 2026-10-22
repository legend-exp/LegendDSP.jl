# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    extremestats(signal::AbstractSamples, start::Real, stop::Real)
    extremestats(signal::RDWaveform, start::RealQuantity, stop::RealQuantity)

Get the extrema and their time positions on `signal` in the interval (`start`,`stop`).
Output is in the format `(min = .., max = .., tmin = .., tmax = ..)`

"""
function extremestats end
export extremestats

function extremestats(input::RadiationDetectorDSP.SamplesOrWaveform, start::RadiationDetectorDSP.RealQuantity, stop::RadiationDetectorDSP.RealQuantity)
    X_axis, Y = RadiationDetectorDSP._get_axis_and_signal(input)
    first_x, step_x = first(X_axis), step(X_axis)
    from = round(Int, ustrip(NoUnits, (start - first_x) / step_x)) + firstindex(X_axis)
    until = round(Int, ustrip(NoUnits, (stop - first_x) / step_x)) + firstindex(X_axis)
    _extremestats_impl(X_axis, Y, from:until)
end

extremestats(input::RadiationDetectorSignals.RDWaveform)  = extremestats(input, first(input.time), last(input.time))
extremestats(input::RadiationDetectorDSP.AbstractSamples) = extremestats(input, firstindex(input), lastindex(input))

function _extremestats_impl(X::AbstractArray{<:RealQuantity}, Y::AbstractArray{<:RealQuantity}, idxs::AbstractUnitRange{<:Integer})
    @assert axes(X) == axes(Y)
    @assert firstindex(X) <= first(idxs) <= last(idxs) <= lastindex(X)
    @assert firstindex(Y) <= first(idxs) <= last(idxs) <= lastindex(Y)

    # `minimum`/`maximum` vectorize, `findmin`/`findmax` do not; `findfirst(isequal(...))`
    # returns the same index as `findmin`/`findmax`: the first extremum, or the first NaN.
    Y_window = view(Y, idxs)
    minv, maxv = minimum(Y_window), maximum(Y_window)
    minidx = idxs[something(findfirst(isequal(minv), Y_window))]
    maxidx = idxs[something(findfirst(isequal(maxv), Y_window))]

    (
        min = minv,
        max = maxv,
        tmin = X[minidx],
        tmax = X[maxidx]
    )
end