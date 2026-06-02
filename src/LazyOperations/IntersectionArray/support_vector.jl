"""
    σ(d::AbstractVector, ia::IntersectionArray)

Return a support vector of an intersection of a finite number of sets in a given
direction.

### Input

- `d`  -- direction
- `ia` -- intersection of a finite number of sets

### Output

A support vector in the given direction.
If the direction has norm zero, the result depends on the individual sets.

### Algorithm

This implementation computes the concrete intersection, which can be expensive.
"""
@validate function σ(d::AbstractVector, ia::IntersectionArray)
    X = concretize(ia)
    return σ(d, X)
end
