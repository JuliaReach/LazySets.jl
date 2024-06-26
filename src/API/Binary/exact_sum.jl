"""
    exact_sum(X::LazySet, Y::LazySet)

Compute the exact sum of two parametric sets.

### Input

- `X` -- parametric set
- `Y` -- parametric set

### Output

A set representing the exact sum, sometimes written ``X ⊞ Y``.
"""
function exact_sum(::LazySet, ::LazySet) end
