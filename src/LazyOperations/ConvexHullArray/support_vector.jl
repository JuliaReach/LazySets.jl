"""
    σ(d::AbstractVector, cha::ConvexHullArray)

Return a support vector of a convex hull of a finite number of sets in a given
direction.

### Input

- `d`   -- direction
- `cha` -- convex hull of a finite number of sets

### Output

A support vector in the given direction.
"""
@validate function σ(d::AbstractVector, cha::ConvexHullArray)
    @assert !isempty(cha.array) "an empty convex hull is not allowed"
    return _σ_union(d, array(cha))
end
