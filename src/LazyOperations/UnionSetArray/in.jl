"""
    in(x::AbstractVector, cup::UnionSetArray)

Check whether a given point is contained in the union of a finite number of
sets.

### Input

- `x`   -- point/vector
- `cup` -- union of a finite number of sets

### Output

`true` iff ``x ∈ cup``.
"""
@validate function in(x::AbstractVector, cup::UnionSetArray)
    return any(X -> x ∈ X, array(cup))
end
