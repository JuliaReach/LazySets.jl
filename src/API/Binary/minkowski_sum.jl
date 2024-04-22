"""
    minkowski_sum(X::LazySet, Y::LazySet)

Compute the Minkowski sum of two sets.

### Input

- `X` -- set
- `Y` -- set

### Output

A set representing the Minkowski sum ``X ⊕ Y``.

### Notes

The Minkowski sum of two sets ``X`` and ``Y`` is defined as

```math
    X ⊕ Y = \\{x + y \\mid x ∈ X, y ∈ Y\\}.
```
"""
function minkowski_sum(::LazySet, ::LazySet) end
