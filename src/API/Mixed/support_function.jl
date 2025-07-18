"""
    ρ(d::AbstractVector, X::LazySet)

Evaluate the support function of a set in a given direction.

### Input

- `d` -- direction
- `X` -- set

### Output

The evaluation of the support function of `X` in direction `d`.

### Notes

The convenience alias `support_function` is also available.

We have the following identity based on the support vector ``σ``:

```math
    ρ(d, X) = d ⋅ σ(d, X)
```
"""
function ρ(::AbstractVector, ::LazySet) end

const support_function = ρ
