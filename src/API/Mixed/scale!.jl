"""
    scale!(α::Real, X::LazySet)

Scale a set by modifying it.

### Input

- `α` -- scalar
- `X` -- set

### Output

The scaled set representing ``α ⋅ X``.
"""
function scale!(::Real, ::LazySet) end
