"""
    dim(cha::ConvexHullArray)

Return the dimension of the convex hull of a finite number of sets.

### Input

- `cha` -- convex hull of a finite number of sets

### Output

The ambient dimension of the convex hull of a finite number of sets, or `0` if
there is no set in the array.
"""
function dim(cha::ConvexHullArray)
    return length(cha.array) == 0 ? 0 : dim(cha.array[1])
end
