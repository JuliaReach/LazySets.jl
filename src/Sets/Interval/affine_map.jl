@validate function affine_map(M, X::Interval, v::AbstractVector; kwargs...)
    @inbounds if length(v) == 1
        return Interval(M[1, 1] * X.dat + v[1])
    end
    return translate(linear_map(M, X; kwargs...), v)
end
