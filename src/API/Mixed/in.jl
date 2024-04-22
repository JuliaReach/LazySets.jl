"""
    ∈(x::AbstractVector, X::LazySet)

Check whether a point lies in a set.

### Input

- `x` -- point/vector
- `X` -- set

### Output

`true` iff ``x ∈ X``.
"""
function ∈(::AbstractVector, ::LazySet) end
