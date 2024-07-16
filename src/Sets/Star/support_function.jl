"""
    ρ(d::AbstractVector, X::Star)

Evaluate the support function of a star.

### Input

- `d` -- direction
- `X` -- star

### Output

The support function in the given direction.
"""
function ρ(d::AbstractVector, X::Star)
    return ρ(At_mul_B(basis(X), d), predicate(X)) + dot(d, center(X))
end
