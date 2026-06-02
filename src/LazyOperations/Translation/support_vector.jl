"""
    σ(d::AbstractVector, tr::Translation)

Return a support vector of a translation.

### Input

- `d`  -- direction
- `tr` -- translation of a set

### Output

A support vector in the given direction.
If the direction has norm zero, the result depends on the wrapped set.
"""
@validate function σ(d::AbstractVector, tr::Translation)
    return tr.v + σ(d, tr.X)
end
