@validate_commutative function distance(x::AbstractVector, H::HalfSpace; p::Real=2)
    @assert length(x) == dim(H) "incompatible dimensions $(length(x)) and $(dim(H))"

    if p != 2
        throw(ArgumentError("`distance` is only implemented for Euclidean norm"))
    end

    N = promote_type(eltype(x), eltype(H))
    a, b = _normalize_halfspace(H, N(2))
    return max(dot(x, a) - b, zero(N))
end
