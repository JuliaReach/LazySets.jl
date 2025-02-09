"""
# Extended help

    ρ(d::AbstractVector, U::Universe)

### Algorithm

If the direction is all zero, the result is zero.
Otherwise, the result is `Inf`.
"""
function ρ(d::AbstractVector, U::Universe)
    @assert length(d) == dim(U) "incompatible dimensions $(length(d)) and $(dim(U))"

    N = promote_type(eltype(d), eltype(U))
    return iszero(d) ? zero(N) : N(Inf)
end
