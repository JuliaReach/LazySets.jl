@commutative function distance(x::AbstractVector, H::HalfSpace)
    N = promote_type(eltype(x), eltype(H))
    a, b = _normalize_halfspace(H, N(2))
    return max(dot(x, a) - b, zero(N))
end
