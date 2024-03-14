# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

struct HaarWaveletFilter <: AbstractRadSigFilter{LinearFiltering}
    Δ::Int
end

export HaarWaveletFilter

struct HaarWaveletFilterInstance{T} <: AbstractRadSigFilterInstance{LinearFiltering}
    Δ::Int
    n_input::Int
end

function fltinstance(flt::HaarWaveletFilter, input::SamplingInfo{T}) where T
    HaarWaveletFilterInstance{T}(flt.Δ, length(input.axis))
end

flt_output_length(fi::HaarWaveletFilterInstance) = (fi.n_input - 1)÷fi.Δ
flt_output_smpltype(fi::HaarWaveletFilterInstance) = flt_input_smpltype(fi)
flt_input_smpltype(::HaarWaveletFilterInstance{T}) where T = T


function rdfilt!(output::AbstractVector{T}, fi::HaarWaveletFilterInstance{T}, input::AbstractVector{T}) where T
    invsq2 = inv(sqrt(T(2)))
    @assert flt_output_length(fi) == length(output) "output length not compatible with filter"
    @inbounds @simd for i in eachindex(output)
        _i = i*fi.Δ
        output[i] = (input[_i] + input[_i + 1]) * invsq2
    end
    output
end