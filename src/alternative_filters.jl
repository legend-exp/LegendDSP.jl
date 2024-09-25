"""
    struct WeightedSavitzkyGolayFilter{T<:RealQuantity} <: AbstractRadFIRFilter

A [Weighted-Savitzky-Golay filter](https://doi.org/10.1021/acsmeasuresciau.1c00054).
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

const WEIGHTTYPES = keys(SGW_COEFFS)

function fltinstance(flt::WeightedSavitzkyGolayFilter, fi::SamplingInfo{T}
    ) where {T}

    fltlen = round(Int, ustrip(NoUnits, flt.length / step(fi.axis)))
    m = (fltlen - 1) ÷ 2

    l = length(fi.axis)
    l_min = T(2*m + 1)
    msg1 = "data too short; min length: $(l_min)"
    @assert length(fi.axis) >= 2*m + 1 msg1
    
    d_max = 2*m
    msg2 = "degree too big for given kernel size; max degree: $(d_max)"
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
        sum = zero(T)
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
        sum = zero(T)
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

function weightFunctionScale(missingFrac::Real, weightType::Int)
    coeffs = SGW_COEFFS[weightType]
    missingFrac <= 0 ? 1 : _w(missingFrac, coeffs)
end

function weight_function(weight_type::Int, x::T) where {T <: Real}
    if x <= -0.999999999999 || x >= 0.999999999999
        return zero(T)
    elseif weight_type == 0
        return one(T)
    elseif weight_type == 1
        decay = T(2.0)
        return exp(-x^2 * decay) + exp(-(x - 2)^2 * decay) + exp(-(x + 2)^2 * decay) -
               2 * exp(-decay) - exp(-9 * decay) # Gaussian-like alpha=2
    elseif weight_type == 2
        return cos(π/2 * x)^2 # Hann
    elseif weight_type == 3
        return (cos(π/2 * x))^4 # Hann-squared
    else weight_type == 4
        return (cos(π/2 * x))^6 # Hann-cube
    end
end

@inline _w(x::Real, v::AbstractVector) = 1 - v[1]/(1 + v[2]*x^v[3])


"""
    struct ModifiedSincFilter{T<:RealQuantity} <: AbstractRadFIRFilter{}

A [Modified-Sinc filter](https://doi.org/10.1021/acsmeasuresciau.1c00054)
Working example
```julia
using RadiationDetectorSignals
using Unitful

n = 600
noise = 1.
t = range(0u"μs", 20u"μs", 2*n)
signal = vcat(zeros(n), 10*ones(n)) + (noise*rand(2*n) .- noise/2)
wf = RDWaveform(t, signal)

# define filter parameters and filter
flt = ModifiedSincFilter(d=4, m=1u"μs")

# apply filter to signal
wf_new = flt(wf)
```

Constructors:

* ```$(FUNCTIONNAME)(; fields...)```

Fields:

$(TYPEDFIELDS)
"""
Base.@kwdef struct ModifiedSincFilter{T<:RealQuantity} <: AbstractRadFIRFilter
    "degree of the filter determining the number of extrema in the kernel"
    d::Int = 2
    "half-width of the kernel"
    m::T = 3
    function ModifiedSincFilter(d::Int, m::T) where T
        s1 = "degree of ModifiedSincFilter needs to be even and be in the interval [2, $(D_MAX)]"
        @assert valid_degree(d) s1

        new{T}(d, m)
    end
end

export ModifiedSincFilter

const D_MAX = 10

valid_degree(d::Int) = iseven(d) && (2 <= d <= D_MAX)

struct ModifiedSincFilterInstance{T<:RealQuantity} <: AbstractRadSigFilterInstance{LinearFiltering}
    d::Int
    m::Int 
    n_input::Int
end

Adapt.adapt_structure(to, flt::ModifiedSincFilter) = flt

function fltinstance(flt::ModifiedSincFilter, fi::SamplingInfo{T}) where {T}
    _m = round(Int, ustrip(NoUnits, flt.m/step(fi.axis)))
    m_min = flt.d/2 + 2
    s2 = "size of kernel too small for given degree"
    @assert _m >= m_min s2

    ModifiedSincFilterInstance{T}(flt.d, _m, length(fi.axis))
end

RadiationDetectorDSP.flt_output_smpltype(fi::ModifiedSincFilterInstance) = RadiationDetectorDSP._floattype(flt_input_smpltype(fi))
RadiationDetectorDSP.flt_input_smpltype(::ModifiedSincFilterInstance{T}) where {T} = T
RadiationDetectorDSP.flt_output_length(fi::ModifiedSincFilterInstance) = fi.n_input
RadiationDetectorDSP.flt_output_time_axis(::ModifiedSincFilterInstance, time::AbstractVector{<:RealQuantity}) = time

function RadiationDetectorDSP.rdfilt!(output::AbstractVector, 
    fi::ModifiedSincFilterInstance, input::AbstractVector)

    @assert firstindex(input) == firstindex(output) "output and input signal don't have matching first indices"
    @assert lastindex(input) == lastindex(output) "output and input signal don't have matching last indices"

    T = eltype(output)
    _y = similar(output, fi.n_input + fi.m*2)
    _unsafe_extendData!(_y, input, fi.m, fi.d)
    kernel = _makeKernel(fi.d, fi.m)
    sum::T = zero(T)
    @inbounds for i in eachindex(input)
        center = fi.m + i
        sum = fma(kernel[1], _y[center], sum)
        @simd for j in Base.OneTo(length(kernel)-2)
            k_j = kernel[j+1]
            sum = fma(k_j, _y[center + j+1], sum)
            sum = fma(k_j, _y[center - (j+1)], sum)
        end
        output[i] = sum
        sum = zero(T)
    end
    output
end

const MS_COEFFS = Dict(
    0 => Float64[],
    2 => Float64[],
    4 => Float64[],
    6 => [0.001717576 0.02437382 1.64375],
    8 => [0.0043993373 0.088211164 2.359375; 0.006146815  0.024715371 3.6359375],
    10 => [0.0011840032 0.04219344 2.746875; 0.0036718843 0.12780383 2.7703125]
)

@inbounds function _makeKernel(d::Int, m::Int)
    coeffs = MS_COEFFS[d]
    T = RadiationDetectorDSP._floattype(Int)
    kernel = Vector{T}(undef, m+1)
    sincArg = (d+4)/2
    sum::T = zero(T)
    κ::Vector{T} = map(n -> κ_j(m, view(coeffs, n, 1:3)), axes(coeffs, 1))
    ν = isodd(d/2) ? 1 : 2
    α = 4
    for i in eachindex(kernel)
        x = (i - 1)/(m + 1)
        k_i = sinc(sincArg*x)
        @simd for j in eachindex(κ)
            s = x*sin((2*j + ν)*π*x)
            k_i = fma(κ[j], s, k_i)
        end
        k_i *= _w(x, α)
        kernel[i] = k_i
        sum += k_i
        i > 1 && (sum += k_i)   # offcenter kernel elements appear twice
    end
    kernel ./= sum
    return kernel
end

"""
    makeFitWeights(d::Int, m::Int) where {U}

return the weights for the linear fit user for linear extrapolation at 
at the right boundary.
"""
function makeFitWeights(d::U, m::Int) where {U <: Integer}
    firstZero = (m+1)/(1.5 + d/2)
    β = 0.7 + 0.14*exp(-0.6*(d-4))
    l = ceil(Int, firstZero*β)
    T = RadiationDetectorDSP._floattype(U)
    w = Vector{T}(undef, l)
    a = pi/2/(firstZero*β)
    for i in eachindex(w)
        w[i] = √(cos(a*(i-1)))
    end
    w
end


function _unsafe_extendData!(y::AbstractVector, x::AbstractVector, m::Int, d::Int)
    lenx = length(x)
    copy!(view(y, m+1:lenx+m), x)
    weights = makeFitWeights(d, m)
    L = min(lenx, length(weights))
    left_x = view(x, 1:L)
    right_x = view(x, lenx-L+1:lenx)
    b1, m1 = weighted_linear_reg(view(weights, 1:L), left_x)
    b2, m2 = weighted_linear_reg(view(weights, 1:L), right_x)
    for p in 1:m
        y[m-p+1] = b1 + m1*(-p)
        y[lenx+m+p] = b2 + m2*(p+L-1)
    end
    y
end


@inline _w(x, α) = exp(-α*x^2) + exp(-α*(x + 2)^2) + exp(-α*(x - 2)^2) - 2*exp(-α) - exp(-9*α)

@inline κ_j(m::Int, v_n::AbstractArray) = v_n[1] + v_n[2]/(v_n[3] - m)^3

"""
    weighted_linear_reg(w::AbstractVector, y::AbstractVector)

Do a weighted linear regression of the data `y` with weights `w` and return
the offset and slope. its assumed that x is a range from 1 to `length(y)`
"""
function weighted_linear_reg(w::AbstractVector, y::AbstractVector)
    T = RadiationDetectorDSP._floattype(promote_type(eltype(w), eltype(y)))
    z = zero(T)

    sum_w = z
    sum_X = z
    sum_Y = z
    sum_X2 = z
    sum_Y2 = z
    sum_XY = z
    @inbounds @fastmath @simd for i in eachindex(y)
        y_i, w_i = y[i], w[i]
        sum_w += w_i
        sum_X = fma(w_i, i-1, sum_X)
        sum_Y = fma(w_i, y_i, sum_Y)
        sum_X2 = fma(w_i, fma(i-1, i-1, 0), sum_X2)
        sum_Y2 = fma(w_i, fma(y_i, y_i, 0), sum_Y2)
        sum_XY = fma(w_i, fma(y_i, i-1, 0), sum_XY)
    end

    inv_sum_w = 1/sum_w
    var_X = sum_X2 - sum_X*sum_X*inv_sum_w
    slope = if sum_w > 0
        s = (sum_XY - sum_X*sum_Y*inv_sum_w)/var_X
        ifelse(isnan(s), zero(T), s)
    else
        T(NaN)
    end
    offset = (sum_Y - slope*sum_X)*inv_sum_w
    offset, slope
end

"""

    struct WhittakerHendersonFilter <: AbstractRadSigFilter{LinearFiltering}

A [Whittaker-Henderson filter](https://doi.org/10.1021/acsmeasuresciau.1c00054).
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
    banded_cholesky_solve!(output, A, input)
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
    banded_cholesky_solve!(x::AbstractVector, A::AbstractMatrix, y::AbstractVector)

solve the linear equation `AA'x = y` and store the result in `x`, where
`A` is the cholesky decomposition of a banded centro symmetric matrix.
`A[1, i]` is the diagonal, `A[2, i]` is the first subdiagonal and so on.
"""
function banded_cholesky_solve!(x::AbstractVector, A::AbstractMatrix, y::AbstractVector)
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
        sum = zero(T)
        @simd for j in (i+1):min(i+dmax, n)
            sum = fma(A[j-i+1, i], x[j], sum)
        end
        x[i] = (x[i] - sum) / A[1, i]
    end
    x
end