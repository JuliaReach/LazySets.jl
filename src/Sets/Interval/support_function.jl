@validate function ρ(d::AbstractVector, X::Interval)
    N = promote_type(eltype(d), eltype(X))
    return @inbounds d[1] * (d[1] > zero(N) ? _max(X) : _min(X))
end
