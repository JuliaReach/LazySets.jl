@commutative function distance(x::AbstractVector, H::Hyperplane)
    @assert length(x) == dim(H) "a vector of length $(length(x)) is " *
                                "incompatible with a set of dimension $(dim(H))"

    N = promote_type(eltype(x), eltype(H))
    a, b = _normalize_halfspace(H, N(2))
    return abs(dot(x, a) - b)
end
