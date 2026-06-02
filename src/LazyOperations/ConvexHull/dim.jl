"""
    dim(ch::ConvexHull)

Return the dimension of a convex hull of two sets.

### Input

- `ch` -- convex hull of two sets

### Output

The ambient dimension of the convex hull of two sets.
"""
function dim(ch::ConvexHull)
    return dim(ch.X)
end
