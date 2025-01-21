@commutative function distance(x::AbstractVector, H::HalfSpace; p::Real=2)
    @assert length(x) == dim(H) "incompatible dimensions $(length(x)) and $(dim(H))"

    N = promote_type(eltype(x), eltype(H))
    if x ∈ H
        return zero(N)
    end
    # formula for hyperplane (wrong!)
    return norm(abs(dot(x, H.a) - H.b) / dot(H.a, H.a) * H.a, p)
end
