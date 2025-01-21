@commutative function distance(x::AbstractVector, H::HalfSpace)
    @assert length(x) == dim(H) "a vector of length $(length(x)) is " *
                                "incompatible with a set of dimension $(dim(H))"

    N = promote_type(eltype(x), eltype(H))
    a, b = _normalize_halfspace(H, N(2))
    return max(dot(x, a) - b, zero(N))
end
