"""
# Extended help

    ρ(d::AbstractVector, U::Universe)

### Algorithm

If the direction is all zero, the result is zero.
Otherwise, the result is `Inf`.
"""
@validate function ρ(d::AbstractVector, U::Universe)
    N = promote_type(eltype(d), eltype(U))
    return iszero(d) ? zero(N) : N(Inf)
end
