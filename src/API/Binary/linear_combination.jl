"""
    linear_combination(X::LazySet, Y::LazySet)

Compute the linear combination of two sets.

### Input

- `X` -- set
- `Y` -- set

### Output

A set representing the linear combination of `X` and `Y`.

### Notes

The linear combination of two sets ``X`` and ``Y`` is defined as

```math
    \\left\\{\\frac{1}{2}(1+λ)x + \\frac{1}{2}(1-λ)y \\mid x ∈ X, y ∈ Y, λ ∈ [-1, 1]\\right\\}.
```
"""
function linear_combination(::LazySet, ::LazySet) end
