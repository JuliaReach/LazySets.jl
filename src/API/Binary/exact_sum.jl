"""
    exact_sum(X::LazySet, Y::LazySet)

Compute the exact sum of two parametric sets.

### Input

- `X` -- parametric set
- `Y` -- parametric set

### Output

A set representing the exact sum ``X ⊞ Y``.

### Notes

The convenience alias `⊞` is also available, which can be typed by
`\\boxplus<tab>`.
"""
function exact_sum(::LazySet, ::LazySet) end

"""
    ⊞(X::LazySet, Y::LazySet)

Convenience alias for the (concrete) `exact_sum` function.

### Notes

"`⊞`" can be typed by `\\boxplus<tab>`.
"""
const ⊞ = exact_sum
