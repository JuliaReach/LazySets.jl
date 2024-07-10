"""
    ρ(d::AbstractVector, U::Universe)

Return the support function of a universe.

### Input

- `d` -- direction
- `U` -- universe

### Output

The support function in the given direction.

### Algorithm

If the direction is all zero, the result is zero.
Otherwise, the result is `Inf`.
"""
function ρ(d::AbstractVector, U::Universe)
    N = promote_type(eltype(d), eltype(U))
    return iszero(d) ? zero(N) : N(Inf)
end
