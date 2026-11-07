# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    extremestats(signal::SamplesOrWaveform, start::RealQuantity, stop::RealQuantity)
    extremestats(signal::RDWaveform)
    extremestats(signal::AbstractSamples)

Minimum, maximum and their time positions of `signal` between `start` and `stop`.

The interval includes the samples nearest to `start` and `stop`. The methods without an interval use the whole signal.

# Returns
- named tuple `(min, max, tmin, tmax)`
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
    
    minv, minidx = findmin(view(Y, idxs))
    maxv, maxidx = findmax(view(Y, idxs))
    
    (
        min = minv,
        max = maxv,
        tmin = X[minidx + first(idxs) - begin],
        tmax = X[maxidx + first(idxs) - begin]
    )

end