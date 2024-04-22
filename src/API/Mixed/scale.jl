
"""
    scale(α::Real, X::LazySet)

Compute the scaling of a set.

### Input

- `α` -- scalar
- `X` -- set

### Output

A set representing ``α ⋅ X``.
"""
function scale(::Real, ::LazySet) end
