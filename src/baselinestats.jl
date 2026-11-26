# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    madstats(signal::SamplesOrWaveform, min::RealQuantity=-Inf, max::RealQuantity=Inf)

Median and width of all samples of `signal` that lie between the lower bound `min` and the upper bound `max`: the width is their median absolute deviation (MAD) from the median, scaled by 1.4826 to estimate the standard deviation of normally distributed data. Samples far from the others, such as a pulse on a trace, hardly move either value.

The default bounds include all samples. Both values are zero if no sample lies between the bounds.

# Returns
- named tuple `(median, sigma)`
"""
function madstats end
export madstats

function madstats(input::RadiationDetectorDSP.SamplesOrWaveform{T}, min::RadiationDetectorDSP.RealQuantity = -Inf * unit(T), max::RadiationDetectorDSP.RealQuantity = Inf * unit(T)) where T
    _, Y = RadiationDetectorDSP._get_axis_and_signal(input)
    _madstats_impl(Y, T(min), T(max))
end

function _madstats_impl(Y::AbstractArray{T}, _min::T, _max::T) where {T <: RealQuantity}
    Y_filt = [y for y in Y if _min <= y <= _max]
    isempty(Y_filt) && return (median = zero(float(T)), sigma = zero(float(T)))
    med = median(Y_filt)
    @inbounds @fastmath for i in eachindex(Y_filt)
        Y_filt[i] = abs(Y_filt[i] - med)
    end
    # 1/Φ⁻¹(3/4) ≈ 1.4826, consistency constant for MAD under normality (Rousseeuw & Croux, J. Amer. Statist. Assoc. 88, 1993)
    (median = float(med), sigma = oftype(float(one(T)), 1.4826) * median!(Y_filt))
end


"""
    sigmaclipstats(signal::SamplesOrWaveform, n_σ::Real)

Mean and width of the samples of `signal` after sigma clipping: the width is the median absolute deviation (MAD) of all samples from their median, scaled by 1.4826 to estimate the standard deviation of normally distributed data, and the mean is taken over the samples within `n_σ` widths of the median.

The mean resolves a level below one ADC count, which the median of the integer samples of a raw waveform does not, and the cut around the median keeps single outliers out of it. `signal` is the stretch itself, for example the samples before a pulse.

# Returns
- named tuple `(mean, sigma)`
"""
function sigmaclipstats end
export sigmaclipstats

function sigmaclipstats(input::RadiationDetectorDSP.SamplesOrWaveform, n_σ::Real)
    _, Y = RadiationDetectorDSP._get_axis_and_signal(input)
    med = median(Y)
    sigma = oftype(float(one(eltype(Y))), 1.4826) * median(abs(y - med) for y in Y)
    (mean = mean(y for y in Y if abs(y - med) <= n_σ * sigma), sigma = sigma)
end
