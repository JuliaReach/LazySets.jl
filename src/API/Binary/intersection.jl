"""
    intersection(X::LazySet, Y::LazySet)

Compute the intersection of two sets.

### Input

- `X` -- set
- `Y` -- set

### Output

A set representing the intersection ``X ∩ Y``.

### Notes

The intersection of two sets ``X`` and ``Y`` is defined as

```math
    X ∩ Y = \\{x \\mid x ∈ X \\text{ and } x ∈ Y\\}.
```
"""
function intersection(::LazySet, ::LazySet) end
