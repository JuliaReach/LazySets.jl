"""
    isbounded(cha::ConvexHullArray)

Check whether a convex hull of a finite number of sets is bounded.

### Input

- `cha` -- convex hull of a finite number of sets

### Output

`true` iff all wrapped sets are bounded.
"""
function isbounded(cha::ConvexHullArray)
    isboundedtype(typeof(cha)) && return true

    return all(isbounded, cha.array)
end
