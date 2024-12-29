@commutative function distance(x::AbstractVector, H::Hyperplane; p::Real=2)
    if p != 2
        throw(ArgumentError("`distance` is only implemented for Euclidean norm"))
    end

    N = promote_type(eltype(x), eltype(H))
    a, b = _normalize_halfspace(H, N(2))
    return abs(dot(x, a) - b)
end
