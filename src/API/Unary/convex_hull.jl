"""
    convex_hull(X::LazySet)

Compute the convex hull of a set.

### Input

- `X` -- set

### Output

A set representing the convex hull of `X`.

### Notes

The convex hull of a set ``X`` is defined as

```math
    \\{λx + (1-λ)y \\mid x, y ∈ X, λ ∈ [0, 1]\\}.
```
"""
function convex_hull(::LazySet) end
