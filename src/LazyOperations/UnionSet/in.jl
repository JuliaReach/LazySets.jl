"""
    in(x::AbstractVector, cup::UnionSet)

Check whether a given point is contained in the union of two sets.

### Input

- `x`   -- point/vector
- `cup` -- union of two sets

### Output

`true` iff ``x ∈ cup``.
"""
@validate function in(x::AbstractVector, cup::UnionSet)
    return x ∈ cup.X || x ∈ cup.Y
end
