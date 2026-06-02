"""
    ρ(d::AbstractVector, ch::ConvexHull)

Evaluate the support function of the convex hull of two sets in a given
direction.

### Input

- `d`  -- direction
- `ch` -- convex hull of two sets

### Output

The evaluation of the support function of the convex hull in the given
direction.
"""
@validate function ρ(d::AbstractVector, ch::ConvexHull)
    return max(ρ(d, ch.X), ρ(d, ch.Y))
end
