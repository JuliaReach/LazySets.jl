"""
# Extended help

    σ(d::AbstractVector, U::Universe)

### Output

A vector with infinity values, except in dimensions where the direction is zero.
"""
function σ(d::AbstractVector, U::Universe)
    @assert length(d) == dim(U) "incompatible dimensions $(length(d)) and $(dim(U))"

    N = promote_type(eltype(d), eltype(U))
    return [iszero(v) ? v : v > zero(N) ? N(Inf) : N(-Inf) for v in d]
end
