# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

struct HaarTypeWaveletFilter <: AbstractRadSigFilter{LinearFiltering}
    down_sampling_rate::Int
    length::Int
    HaarTypeWaveletFilter(down_sampling_rate::Int) = new(down_sampling_rate, 2)
end

export HaarTypeWaveletFilter

struct HaarTypeWaveletFilterInstance{T} <: AbstractRadSigFilterInstance{LinearFiltering}
    down_sampling_rate::Int
    length::Int
    n_input::Int
end

function RadiationDetectorDSP.fltinstance(flt::HaarTypeWaveletFilter, input::SamplingInfo{T}) where T
    HaarTypeWaveletFilterInstance{T}(flt.down_sampling_rate, flt.length, length(input.axis))
end

RadiationDetectorDSP.flt_output_smpltype(fi::HaarTypeWaveletFilterInstance) = flt_input_smpltype(fi)
RadiationDetectorDSP.flt_input_smpltype(::HaarTypeWaveletFilterInstance{T}) where T = T
RadiationDetectorDSP.flt_output_length(fi::HaarTypeWaveletFilterInstance) = ceil(Int, fi.n_input/fi.down_sampling_rate)
RadiationDetectorDSP.flt_output_time_axis(fi::HaarTypeWaveletFilterInstance, time::AbstractVector{<:RealQuantity}) = time[1:fi.down_sampling_rate:end]

function RadiationDetectorDSP.rdfilt!(output::AbstractVector{T}, fi::HaarTypeWaveletFilterInstance{T}, input::AbstractVector{T}) where T
    invsqrt = inv(sqrt(T(fi.length)))
    @assert flt_output_length(fi) == length(output) "output length not compatible with filter"
    Δidx = fi.length - 1
    output[1] = input[1] / invsqrt
    output[end] = input[end] / invsqrt
    for i in 2:length(output)-1
        j = (i-1)*fi.down_sampling_rate + 1
        output[i] = sum(view(input, j:j+Δidx)) * invsqrt
    end
    output
end