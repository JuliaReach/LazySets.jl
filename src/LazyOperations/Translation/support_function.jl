"""
    ρ(d::AbstractVector, tr::Translation)

Evaluate the support function of a translation.

### Input

- `d`  -- direction
- `tr` -- translation of a set

### Output

The evaluation of the support function in the given direction.
"""
@validate function ρ(d::AbstractVector, tr::Translation)
    return dot(d, tr.v) + ρ(d, tr.X)
end
