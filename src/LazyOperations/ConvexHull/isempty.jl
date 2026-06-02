"""
    isempty(ch::ConvexHull)

Check whether the convex hull of two sets is empty.

### Input

- `ch` -- convex hull

### Output

`true` iff both wrapped sets are empty.
"""
function isempty(ch::ConvexHull)
    return isempty(ch.X) && isempty(ch.Y)
end
