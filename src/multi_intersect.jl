
"""
struct MultiIntersect <: Function
    Finds the x values at which a signal exceeds the provided 
    `threshold_ratios` the first time. 
    Constructors:
    * ```$(FUNCTIONNAME)(; fields...)```
    Fields:
    $(TYPEDFIELDS)
"""
Base.@kwdef struct MultiIntersect{T<:RadiationDetectorDSP.RealQuantity} <: Function 
    "ratios which determine the thresholds"
    threshold_ratios::Vector{Float64} = collect(0.01:0.01:0.9)
    "minimum time-over-threshold"
    mintot::T = 4
    "half window length of polynomial fit"
    n::Int = 1
    "degree of polynomial"
    d::Int = 1
    "upsampling rate"
    sampling_rate::Int = 1
end

export MultiIntersect

function (f::MultiIntersect)(input::RadiationDetectorDSP.SamplesOrWaveform)

    X_axis, Y = RadiationDetectorDSP._get_axis_and_signal(input)
    tmp = f.mintot / step(X_axis)
    min_n_over_thresh = max(1, round(Int, ustrip(NoUnits, tmp)))
    thresholds = f.threshold_ratios .* maximum(Y)
    _find_intersects(
        X_axis, Y, min_n_over_thresh, thresholds, f.n, f.d, f.sampling_rate)
end

function _find_intersects(X::AbstractVector{T}, 
    Y::AbstractVector{U}, min_n_over_thresh::Int, 
    thresholds::Vector{<:RealQuantity}, n::Int, degree::Int, rate::Int
    ) where {T <: RealQuantity, U <: RealQuantity}
    
    l = length(thresholds)
    @assert axes(X) == axes(Y)

    # ToDo: What if eltype(Y) is based on ForwardDiff.Dual, but eltype(X) is not?
    R = float(T)
    intersect_x::Array{R} = zeros(R, l)


    isempty(Y) && return intersect_x

    # GPU-friendly branch-free code:

    cand_pos::Int = firstindex(Y) + 1
    intersect_pos::Array{Int} = fill(firstindex(Y) + 1, l)
    y_high_counter::Int = ifelse(first(Y) >= first(thresholds), min_n_over_thresh + 1, 0)
    intersect_counter::Int = 1
    i::Int = firstindex(Y)
    while (i <= lastindex(Y) && intersect_counter < l + 1)
        y = Y[i]
        y_is_high = y >= thresholds[intersect_counter]
        first_high_y = y_high_counter == 0
        cand_pos = ifelse(y_is_high && first_high_y, i, cand_pos)
        y_high_counter = ifelse(y_is_high, y_high_counter + 1, 0)
        new_intersect_found = y_high_counter == min_n_over_thresh
        pos = ifelse(new_intersect_found, cand_pos, intersect_pos[intersect_counter])
        intersect_pos[intersect_counter] = pos

        # reset running index i to previous found intersect pos
        i = ifelse(new_intersect_found, pos, i + 1)
        intersect_counter = ifelse(new_intersect_found, intersect_counter + 1, intersect_counter)
        y_high_counter = ifelse(new_intersect_found, 0, y_high_counter)
    end
    
    # check lowest and highest boundary to ensure @inbounds for polynomial fit
    outside_left_boundary = !(intersect_pos[begin] - n >= firstindex(Y))
    outside_right_boundary = !(intersect_pos[end] + n - 1 <= lastindex(Y))
    errmsg = "cannot interpolate intersect on left boundary"
    @assert !(outside_left_boundary || outside_right_boundary) errmsg

    A = RadiationDetectorDSP._lsq_fit_matrix(0:2*n - 1, degree)
    m = 2*n*rate
    x_up = range(0, 2*n - 1, m)
    # precompute vandermonde matrix for upsampled x axis
    x_van = x_up .^ (0:degree)'
    V = promote_type(eltype(A), float(U))
    y_up::Vector{V} = zeros(V, m)

    @inbounds for j in eachindex(intersect_pos)
        from = intersect_pos[j] - n
        to = intersect_pos[j] + n - 1
        _x_axis = range(X[from], X[to], m)

        # fit polynomial and save evaluation of it at upsampled position
        _lsqfitatwindow!(A, view(Y, from:to), x_van, y_up)

        # search for threshold in upsampled signal window
        intersect_x[j] = RadiationDetectorDSP._find_intersect_impl(
            _x_axis, y_up, thresholds[j], 1).x

        # reset upsampled values
        y_up .= zero(V)
    end
    return intersect_x
end

function _lsqfitatwindow!(A_slqfit::AbstractMatrix{<:Real},
    Y::AbstractVector{<:Real}, x_vandermonde::AbstractMatrix{<:Real}, 
    res::AbstractVector{R}
    ) where {R <: Real}

    @assert axes(A_slqfit, 1) == axes(Y, 1)
    @assert axes(A_slqfit, 2) == axes(x_vandermonde, 2)
    @assert axes(res, 1) == axes(x_vandermonde, 1)

    @inbounds for j in axes(A_slqfit, 2)
        cⱼ₋₁::R = zero(R)
        @simd for i in axes(A_slqfit, 1)
            cⱼ₋₁ = fma(A_slqfit[i, j], Y[i], cⱼ₋₁)
        end
        @inbounds @simd for idx in eachindex(res)
            res[idx] = fma(cⱼ₋₁, x_vandermonde[idx, j], res[idx])
        end
    end
    res
end