@doc """
    translate(X::LazySet, v::AbstractVector)

Compute the translation of a set with a vector.

### Input

- `X` -- set
- `v` -- vector

### Output

A set representing ``X + \\{v\\}``.
"""
translate(::LazySet, ::AbstractVector)
function translate end
