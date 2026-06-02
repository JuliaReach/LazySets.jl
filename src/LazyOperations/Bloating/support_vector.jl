"""
    σ(d::AbstractVector, B::Bloating)

Return the support vector of a bloated set in a given direction.

### Input

- `d` -- direction
- `B` -- bloated set

### Output

The support vector of the bloated set in the given direction.
"""
@validate function σ(d::AbstractVector, B::Bloating)
    @assert !iszero(d) "the support vector in the zero direction is undefined"
    @assert B.ε >= 0 || B.p > 1 "the support vector for negative bloating " *
                                "in the 1-norm is not implemented"

    return σ(d, B.X) +
           sign_cadlag(B.ε) * σ(d, _bloating_ball(abs(B.ε), B.p, dim(B)))
end
