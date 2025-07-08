"""
    σ(d::AbstractVector, X::LazySet)

Compute a support vector of a set in a given direction.

### Input

- `d` -- direction
- `X` -- set

### Output

A support vector of `X` in direction `d`.

### Notes

The convenience alias `support_vector` is also available.
"""
function σ(::AbstractVector, ::LazySet) end

const support_vector = σ
