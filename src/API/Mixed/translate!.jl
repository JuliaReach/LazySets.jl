"""
    translate!(X::LazySet, v::AbstractVector)

Translate a set with a vector by modifying it.

### Input

- `X` -- set
- `v` -- vector

### Output

The translated set representing ``X + \\{v\\}``.
"""
function translate!(::LazySet, ::AbstractVector) end
