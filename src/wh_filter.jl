"""

    WhittakerHendersonFilter <: AbstractRadSigFilter{LinearFiltering}

A [Whittaker-Henderson filter](https://doi.org/10.1021/acsmeasuresciau.3c00017).
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
flt = WhittakerHendersonFilter(p=3, λ=4)

# apply filter to signal
wf_new = flt(wf)
"""
Base.@kwdef struct WhittakerHendersonFilter <: AbstractRadSigFilter{LinearFiltering}
    "highest differentiation degree"
    p::Int = 1

    "penalty strengh imposed on non-smooth functions"
    λ::Real = 1
end

export WhittakerHendersonFilter

struct WhittakerHendersonFilterInstance{T<:RealQuantity} <: AbstractRadSigFilterInstance{LinearFiltering}
    p::Int
    λ::Real
    n_input::Int
end

Adapt.adapt_structure(to, flt::WhittakerHendersonFilter) = flt

function fltinstance(flt::WhittakerHendersonFilter, fi::SamplingInfo{T}) where {T}
    WhittakerHendersonFilterInstance{T}(flt.p, flt.λ, length(fi.axis))
end

function rdfilt!(output::AbstractVector, fi::WhittakerHendersonFilterInstance, input::AbstractVector)
    A = makeIPlusLambdaDprimeD(fi.λ, fi.p, length(output))
    choleskyL!(A)
    solve!(output, A, input)
end

RadiationDetectorDSP.flt_output_smpltype(fi::WhittakerHendersonFilterInstance{T}) where {T} = RadiationDetectorDSP._floattype(flt_input_smpltype(fi))
RadiationDetectorDSP.flt_input_smpltype(::WhittakerHendersonFilterInstance{T}) where {T} = T
RadiationDetectorDSP.flt_output_length(fi::WhittakerHendersonFilterInstance) = RadiationDetectorDSP.flt_input_length(fi)
RadiationDetectorDSP.flt_input_length(fi::WhittakerHendersonFilterInstance) = fi.n_input
RadiationDetectorDSP.flt_output_time_axis(::WhittakerHendersonFilterInstance, time::AbstractVector) = time

_coeffs(p) = map(n -> binomial(p, n)*(-1)^(n + p), 0:p)
@inline _coeff(p, n) = binomial(p, n-1)*(-1)^(n-1+p)

"""
    makeIPlusLambdaDprimeD(λ::T, p::Int, N::Int) where {T}

build the centro symmetric banded matrix which needs to be inverted 
later: I + λ*D'D
D is the finite difference matrix of order `p`.
"""
function makeIPlusLambdaDprimeD(λ::T, p::Int, N::Int) where {T}
    @assert N > p "Order ($p) must be less than number of points ($N)"

    U = RadiationDetectorDSP._floattype(T)
    out = Array{U}(undef, p+1, N)  # Initialize output matrix
    sum::U = zero(U)
    @inbounds for d in 0:p
        len = N - d
        for i in 1:((len+1) ÷ 2)
            sum = 0.0
            from = max(1, i - len + p - d + 1)
            to = min(i, p - d + 1)
            @simd for j in from:to
                sum = fma(_coeff(p, j), _coeff(p, j+d), sum)
            end
            out[d+1, i] = (1 - min(d, 1)) + λ*sum
            out[d+1, len-i+1] = (1 - min(d, 1)) + λ*sum
        end
    end
    out
end

"""

    choleskyL!(b::AbstractMatrix{T})

inplace cholesky decomposition of a banded symmetric matrix, where 
b[1, i] contains the diagonal elements, b[2, i] the elements of the 
first subdiagonal and so on.
"""
function choleskyL!(b::AbstractMatrix{T}) where {T}
    dmax = size(b, 1) - 1
    @inbounds for i in axes(b, 2)
        for j in max(1, i-dmax):i
            sum = 0.
            for k in max(1, i-dmax):(j-1)
                sum += b[i-k+1, k] * b[j-k+1, k]
            end
            if (i == j)
                sqrtArg = b[1, i] - sum
                @assert sqrtArg > 0 "Matrix is not positive definite"
                b[1, i] = sqrt(sqrtArg)
            else
                b[i-j+1, j]  = 1.0 / b[1, j] * (b[i-j+1, j] - sum)
            end
        end
    end
    b
end

"""
    solve!(x::AbstractVector, A::AbstractMatrix, y::AbstractVector)

solve the lienar equation AA'x = y and store the result in x, where
`A` is the cholesky decomposition of a banded centro symmetric matrix.
A[1, i] is the diagonal, A[2, i] is the first subdiagonal and so on.
"""
function solve!(x::AbstractVector, A::AbstractMatrix, y::AbstractVector)
    n = size(A, 2)
    dmax = size(A, 1) - 1
    T = RadiationDetectorDSP._floattype(eltype(x))
    sum::T = zero(T)

    @inbounds for i in axes(A, 2)
        sum = zero(T)
        @simd for j in max(1, i-dmax):(i-1)
            sum = fma(A[i-j+1, j], x[j], sum)
        end
        x[i] = (y[i] - sum) / A[1, i] 
    end

    # Backward substitution
    @inbounds for i in reverse(axes(A, 2))
        sum = 0.0
        @simd for j in (i+1):min(i+dmax, n)
            sum = fma(A[j-i+1, i], x[j], sum)
        end
        x[i] = (x[i] - sum) / A[1, i]
    end
    x
end
