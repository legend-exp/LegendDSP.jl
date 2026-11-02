# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    struct CenteredConvolutionFilter <: AbstractRadFIRFilter

Linear convolution of the signal with a two-sided kernel `coeffs` of odd length whose middle tap is lag 0, e.g. the
inverse kernel of a SiPM channel. The output has the length and the time axis of the input; the signal is taken as
zero beyond its ends. FFT based, so long kernels cost no more than short ones.

Constructors:

* ```$(FUNCTIONNAME)(; fields...)```
* ```$(FUNCTIONNAME)(coeffs)```

Fields:

$(TYPEDFIELDS)
"""
Base.@kwdef struct CenteredConvolutionFilter{TV <: AbstractVector{<:Real}} <: AbstractRadFIRFilter
    "filter taps, odd length, lag 0 in the middle"
    coeffs::TV
end

export CenteredConvolutionFilter

struct CenteredConvolutionFilterInstance{T} <: AbstractRadSigFilterInstance{LinearFiltering}
    "spectrum of the taps, zero padded to the linear-convolution length"
    rfft_h::Vector{ComplexF64}
    n_taps::Int
    n_input::Int
end

function RadiationDetectorDSP.fltinstance(flt::CenteredConvolutionFilter, si::SamplingInfo{T}) where {T}
    m = length(flt.coeffs)
    isodd(m) || throw(ArgumentError("CenteredConvolutionFilter needs an odd number of taps, got $m"))
    n = _smpllen(si)
    CenteredConvolutionFilterInstance{T}(rfft(vcat(Float64.(flt.coeffs), zeros(n - 1))), m, n)
end

RadiationDetectorDSP.flt_output_smpltype(::CenteredConvolutionFilterInstance{T}) where {T} = _floattype(T)
RadiationDetectorDSP.flt_input_smpltype(::CenteredConvolutionFilterInstance{T}) where {T} = T
RadiationDetectorDSP.flt_output_length(fi::CenteredConvolutionFilterInstance) = fi.n_input
RadiationDetectorDSP.flt_input_length(fi::CenteredConvolutionFilterInstance) = fi.n_input
RadiationDetectorDSP.flt_output_time_axis(::CenteredConvolutionFilterInstance, time::AbstractVector) = time

function rdfilt!(y::AbstractVector, fi::CenteredConvolutionFilterInstance, x::AbstractVector)
    n = fi.n_input
    m = fi.n_taps
    c = irfft(rfft(vcat(Float64.(x), zeros(m - 1))) .* fi.rfft_h, n + m - 1)
    mid = (m + 1) ÷ 2
    y .= @view c[mid:mid+n-1]
    y
end
