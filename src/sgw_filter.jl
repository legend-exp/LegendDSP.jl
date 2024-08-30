"""
    struct WeightedSavitzkyGolayFilter{T<:RealQuantity} <: AbstractRadFIRFilter

A [Weighted-Savitzky-Golay filter](https://doi.org/10.1021/acsmeasuresciau.3c00017).
Working example:
```julia
using RadiationDetectorSignals
using Unitful

n = 600
noise = 1.
t = range(0u"μs", 20u"μs", 2*n)
signal = vcat(zeros(n), 10*ones(n)) + (noise*rand(2*n) .- noise/2)
wf = RDWaveform(t, signal)

# define filter parameters and filter
flt = WeightedSavitzkyGolayFilter(length=1u"μs", degree=3, weightType=2)

# apply filter to signal
wf_new = flt(wf)
```


Constructors:

* ```$(FUNCTIONNAME)(; fields...)```

Fields:

$(TYPEDFIELDS)
"""
Base.@kwdef struct WeightedSavitzkyGolayFilter{T<:RealQuantity} <: AbstractRadFIRFilter
    "total filter length"
    length::T

    "Polynomial degree"
    degree::Int = 1

    "weight function to use"
    weightType::Int = 0

    function WeightedSavitzkyGolayFilter(length::T, degree::Int, 
        weightType::Int) where {T}

        @assert weightType ∈ WEIGHTTYPES "unvalid weight function provided"

        new{T}(length, degree, weightType)
    end
end

export WeightedSavitzkyGolayFilter

struct WeightedSavitzkyGolayFilterInstance{T<:RealQuantity} <: AbstractRadSigFilterInstance{NonlinearFiltering}
    m::Int
    degree::Int
    weightType::Int
    n_input::Int
end

const SGW_COEFFS = Dict(
    0 => [1., 1., -1.],                  # normal Sawitzky-Golay filter
    1 => [0.68096, 0.36358, -3.68528],   # GAUSS2
    2 => [0.67574, 0.35440, -3.61580],   # HANN
    3 => [0.63944, 0.28417, -5.508],     # HANNSQR
    4 => [0.62303, 0.25310, -7.07317]    # HANNCUBE
)

const WEIGHTTYPES = (0, 1, 2, 3, 4)

function fltinstance(flt::WeightedSavitzkyGolayFilter, fi::SamplingInfo{T}
    ) where {T}

    fltlen = round(Int, ustrip(NoUnits, flt.length / step(fi.axis)))
    m = (fltlen - 1) ÷ 2

    l = length(fi.axis)
    l_min = T(2*m + 1)
    msg1 = "data too short; min length: $(l_min)"
    @assert length(fi.axis) >= 2*m + 1 msg1
    
    d_max = 2*m
    msg2 = "degree to big for given kernel size; max degree: $(d_max)"
    @assert flt.degree <= 2*m msg2

    WeightedSavitzkyGolayFilterInstance{T}(m, flt.degree, flt.weightType, l)
end

RadiationDetectorDSP.flt_output_smpltype(fi::WeightedSavitzkyGolayFilterInstance) = RadiationDetectorDSP._floattype(flt_input_smpltype(fi))
RadiationDetectorDSP.flt_input_smpltype(::WeightedSavitzkyGolayFilterInstance{T}) where T = T
RadiationDetectorDSP.flt_output_length(fi::WeightedSavitzkyGolayFilterInstance) = fi.n_input
RadiationDetectorDSP.flt_output_time_axis(::WeightedSavitzkyGolayFilterInstance, time::AbstractVector{<:RealQuantity}) = time

function RadiationDetectorDSP.rdfilt!(output::AbstractVector, 
    fi::WeightedSavitzkyGolayFilterInstance, input::AbstractVector)

    T = eltype(output)
    k_len = fi.m * 2 + 1
    Polynomials = Array{T}(undef, fi.degree+1, k_len)
    Kernel = Array{T}(undef, k_len)
    weights = Vector{T}(undef, k_len)
    _unsafe_convolve!(output, Polynomials, Kernel, weights, input, fi.m, 
        fi.degree, fi.weightType)
end

Adapt.adapt_structure(to, flt::WeightedSavitzkyGolayFilter) = flt

@inbounds function _unsafe_convolve!(y::AbstractVector, P::AbstractMatrix, 
    kernel::AbstractVector, weights::AbstractVector, x::AbstractVector,
    m::Int, d::Int, weightType::Int)

    T = eltype(y)
    sum::T = zero(T)
    pRight::Int = 0
    scale::T = zero(T)
    L = length(x)
    # left near-boundary and interior points
    for i in Base.OneTo(L - m)
        pLeft = min(m, i-1)
        scale = weightFunctionScale((m - pLeft)/m, weightType)
        pRight = floor(Int, (m+1)/scale)
        pRight = ifelse((pRight + pLeft) > 2*m, 2*m - pLeft, pRight)
        len_k = pRight + pLeft + 1
        _unsafe_make_left_kernel!(kernel, P, weights, pLeft, len_k, 
            scale, m, d, weightType)
        for j in Base.OneTo(len_k)
            idx = max(i-m-1, 0) + j
            sum = fma(kernel[j], x[idx], sum)
        end
        y[i] = sum
        sum = 0.
    end
    # near boundary points at the right
    sum = zero(T)
    for i in Base.OneTo(m)
        pLeft = m-i+1
        scale = weightFunctionScale((m - pLeft)/m, weightType)
        pRight = floor(Int, (m+1)/scale)
        pRight = ifelse((pRight + pLeft) > 2*m, 2*m - pLeft, pRight)
        len_k = pRight + pLeft + 1
        _unsafe_make_left_kernel!(kernel, P, weights, pLeft, len_k, 
            scale, m, d, weightType)
        len_k = pRight + pLeft + 1
        for j in Base.OneTo(len_k)
            idx = L + 1 - j
            sum = fma(kernel[j], x[idx], sum)
        end
        y[L-m+i] = sum
        sum = 0.
    end
    y
end

@inbounds function _unsafe_make_left_kernel!(kernel::AbstractVector, 
    polynomials::AbstractMatrix, weights::AbstractVector, pLeft::Int, 
    k_len::Int, scale::T, m::Int, d::Int, weightType::Int) where T

    for i in Base.OneTo(k_len - pLeft)
        w_i = weight_function(weightType, (i-1)*scale/(m+1))
        weights[pLeft+i] = w_i
        ((i != 1) && (i < pLeft)) && (weights[pLeft-i] = w_i)
    end
    inv_sumw = 1.0/sqrt(sum(weights))
    @simd for i in axes(polynomials, 2)
        polynomials[1, i] = inv_sumw
        kernel[i] = 0.
    end
    for o in Base.OneTo(d)
        @simd for i in Base.OneTo(k_len)
            polynomials[o+1, i] = polynomials[o, i] * (i-1-pLeft)
        end
    end

    # modified gram-schmidt orthonormalization
    dot_p_o_p_u::T = 0.
    for o in 2:d+1
        p_o = view(polynomials, o, 1:k_len)
        for u in Base.OneTo(o-1)
            p_u = view(polynomials, u, 1:k_len)
            dot_p_o_p_u = -_unsafe_dot_w(p_u, p_o, weights)
            @simd for i in Base.OneTo(k_len)
                p_o[i] = fma(p_u[i], dot_p_o_p_u, p_o[i])
            end
        end
        norm_p_o = _unsafe_dot_w(p_o, p_o, weights)
        p_o ./= sqrt(norm_p_o)
    end
    
    # build kernel according to:
    # aᵢ = ∑ⱼ∑ᵢ wᵢ⋅pⱼ(i)⋅pⱼ(pleft), where i ∈ 0..n, j ∈ -m..m 
    
    for o in Base.OneTo(d+1)
        s = polynomials[o , pLeft+1]
        @simd for i in Base.OneTo(k_len)
            kernel[i] = fma(polynomials[o, i], weights[i]*s, kernel[i])
        end
    end
    kernel
end

function _unsafe_dot_w(x::AbstractVector, y::AbstractVector, w::AbstractVector)
    dotprod = 0.
    @simd for i in eachindex(x)
        dotprod = fma(x[i], fma(y[i], w[i], 0.), dotprod)
    end
    dotprod
end

function weightFunctionScale(missingFrac, weightType::Int)
    coeffs = SGW_COEFFS[weightType]
    missingFrac <= 0 ? 1 : _w(missingFrac, coeffs)
end

function weight_function(weight_type::Int, x::T) where T
    decay = 2.0 # for GAUSS only
    if x <= -0.999999999999 || x >= 0.999999999999
        return 0.0
    elseif weight_type == 0
        return 1.0
    elseif weight_type == 1
        return exp(-x^2 * decay) + exp(-(x - 2)^2 * decay) + exp(-(x + 2)^2 * decay) -
               2 * exp(-decay) - exp(-9 * decay) # Gaussian-like alpha=2
    elseif weight_type == 2
        return cos(0.5 * π * x)^2 # Hann
    elseif weight_type == 3
        return (cos(0.5 * π * x))^4 # Hann-squared
    else weight_type == 4
        return (cos(0.5 * π * x))^6 # Hann-cube
    end
end

@inline _w(x, v::AbstractVector) = 1 - v[1]/(1 + v[2]*x^v[3])