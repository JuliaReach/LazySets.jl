"""
    ρ(d::AbstractVector, B::Bloating)

Return the support function of a bloated set in a given direction.

### Input

- `d` -- direction
- `B` -- bloated set

### Output

The support function of the bloated set in the given direction.
"""
@validate function ρ(d::AbstractVector, B::Bloating)
    @assert !iszero(d) "the support function in the zero direction is undefined"

    return ρ(d, B.X) +
           sign_cadlag(B.ε) * ρ(d, _bloating_ball(abs(B.ε), B.p, dim(B)))
end
