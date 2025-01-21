@commutative function distance(x::AbstractVector, H::Hyperplane)
    N = promote_type(eltype(x), eltype(H))
    a, b = _normalize_halfspace(H, N(2))
    return abs(dot(x, a) - b)
end
