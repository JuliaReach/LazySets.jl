
@doc """
    scale(α::Real, X::LazySet)

Compute the scaling of a set.

### Input

- `α` -- scalar
- `X` -- set

### Output

A set representing ``α ⋅ X``.
"""
scale(::Real, ::LazySet)
function scale end
