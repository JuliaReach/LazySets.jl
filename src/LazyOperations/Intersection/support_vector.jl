"""
    σ(d::AbstractVector, cap::Intersection)

Return a support vector of an intersection of two sets in a given
direction.

### Input

- `d`   -- direction
- `cap` -- intersection of two sets

### Output

A support vector in the given direction.

### Algorithm

We compute the concrete intersection, which may be expensive.
"""
@validate function σ(d::AbstractVector, cap::Intersection)
    X = concretize(cap)
    return σ(d, X)
end
