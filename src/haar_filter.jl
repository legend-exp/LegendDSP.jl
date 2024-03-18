# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

struct HaarTypeWaveletFilter <: AbstractRadSigFilter{LinearFiltering}
    sampling_rate::Int
    length::Int
end

export HaarTypeWaveletFilter

struct HaarTypeWaveletFilterInstance{T} <: AbstractRadSigFilterInstance{LinearFiltering}
    sampling_rate::Int
    length::Int
    n_input::Int
end

function fltinstance(flt::HaarTypeWaveletFilter, input::SamplingInfo{T}) where T
    HaarTypeWaveletFilterInstance{T}(flt.sampling_rate, flt.length, length(input.axis))
end

flt_output_smpltype(fi::HaarTypeWaveletFilterInstance) = flt_input_smpltype(fi)
flt_input_smpltype(::HaarTypeWaveletFilterInstance{T}) where T = T
flt_output_length(fi::HaarTypeWaveletFilterInstance) = 
    ceil(Int, fi.n_input/fi.sampling_rate)
flt_output_time_axis(fi::HaarTypeWaveletFilter, 
    time::AbstractVector{<:RealQuantity}) = time[1:fi.sampling_rate:end]

function rdfilt!(output::AbstractVector{T}, fi::HaarTypeWaveletFilterInstance{T}, input::AbstractVector{T}) where T
    invsqrt = inv(sqrt(T(fi.length)))
    @assert flt_output_length(fi) == length(output) "output length not compatible with filter"
    @inbounds @simd for i in eachindex(output)
        to = min(fi.n_input, (i-1)*fi.sampling_rate + 1)
        from = max(1, to + 1 - fi.length)
        output[i] = sum(view(input, from:to)) * invsqrt
    end
    output
end