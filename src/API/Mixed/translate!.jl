@doc """
    translate!(X::LazySet, v::AbstractVector)

Translate a set with a vector by modifying it.

### Input

- `X` -- set
- `v` -- vector

### Output

The translated set representing ``X + \\{v\\}``.
"""
translate!(::LazySet, ::AbstractVector)
function translate! end
