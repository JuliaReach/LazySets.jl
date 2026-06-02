"""
    isbounded(ch::ConvexHull)

Check whether the convex hull of two sets is bounded.

### Input

- `ch` -- convex hull of two sets

### Output

`true` iff both wrapped sets are bounded.
"""
function isbounded(ch::ConvexHull)
    return isbounded(ch.X) && isbounded(ch.Y)
end
