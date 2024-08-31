
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