@commutative function distance(x::AbstractVector, H::HalfSpace; p::Real=2)
    @assert length(x) == dim(H) "a vector of length $(length(x)) is " *
                                "incompatible with a set of dimension $(dim(H))"

    if p != 2
        throw(ArgumentError("this method is only implemented for the Euclidean norm"))
    end

    N = promote_type(eltype(x), eltype(H))
    a, b = _normalize_halfspace(H, N(2))
    return max(dot(x, a) - b, zero(N))
end
