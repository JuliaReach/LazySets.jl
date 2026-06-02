"""
    in(x::AbstractVector, tr::Translation)

Check whether a given point is contained in the translation of a set.

### Input

- `x`  -- point/vector
- `tr` -- translation of a set

### Output

`true` iff ``x ∈ tr``.

### Algorithm

This implementation relies on the set-membership function for the wrapped set
`tr.X`, since ``x ∈ X ⊕ v`` iff ``x - v ∈ X``.
"""
@validate function in(x::AbstractVector, tr::Translation)
    return x - tr.v ∈ tr.X
end
