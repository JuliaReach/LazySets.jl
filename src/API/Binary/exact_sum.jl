"""
    exact_sum(X::LazySet, Y::LazySet)

Compute the exact sum of two parametric sets.

### Input

- `X` -- parametric set
- `Y` -- parametric set

### Output

A set representing the exact sum, sometimes written ``X ⊞ Y``.

### Notes

For parametric sets, the exact sum behaves like the Minkowski sum, except that
the parameters are shared. Thus, for nonparametric sets, it coincides with the
Minkowski sum.
"""
function exact_sum(::LazySet, ::LazySet) end
