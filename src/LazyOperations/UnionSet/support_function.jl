"""
    ρ(d::AbstractVector, cup::UnionSet)

Evaluate the support function of the union of two sets in a given direction.

### Input

- `d`   -- direction
- `cup` -- union of two sets

### Output

The evaluation of the support function in the given direction.

### Algorithm

The support function of the union of two sets ``X`` and ``Y`` evaluates to the
maximum of the support-function evaluations of ``X`` and ``Y``.
"""
@validate function ρ(d::AbstractVector, cup::UnionSet)
    X, Y = cup.X, cup.Y
    return max(ρ(d, X), ρ(d, Y))
end
