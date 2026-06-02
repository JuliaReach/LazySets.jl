"""
    σ(d::AbstractVector, ch::ConvexHull)

Return a support vector of the convex hull of two sets in a given direction.

### Input

- `d`  -- direction
- `ch` -- convex hull of two sets

### Output

A support vector of the convex hull in the given direction.
"""
@validate function σ(d::AbstractVector, ch::ConvexHull)
    σ1 = σ(d, ch.X)
    σ2 = σ(d, ch.Y)
    ρ1 = dot(d, σ1)
    ρ2 = dot(d, σ2)
    return ρ1 >= ρ2 ? σ1 : σ2
end
