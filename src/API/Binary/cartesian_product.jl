"""
    cartesian_product(X::LazySet, Y::LazySet)

Compute the Cartesian product of two sets.

### Input

- `X` -- set
- `Y` -- set

### Output

A set representing the Cartesian product ``X × Y``.

### Notes

The Cartesian product of two sets ``X`` and ``Y`` is defined as

```math
    X × Y = \\{[x, y] \\mid x ∈ X, y ∈ Y\\}.
```
"""
function cartesian_product(::LazySet, ::LazySet) end
