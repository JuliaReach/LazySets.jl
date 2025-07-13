@validate function σ(d::AbstractVector, X::Interval)
    N = promote_type(eltype(d), eltype(X))
    return @inbounds d[1] > zero(N) ? high(X) : low(X)
end
