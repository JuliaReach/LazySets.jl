"""
    isconvextype(T::Type{<:LazySet})

Check whether `T` is convex just by using type information.

### Input

- `T` -- subtype of `LazySet`

### Output

`true` iff the set type only represents convex sets.

### Notes

Since this operation only acts on types (not on values), it can return false
negatives, i.e., there may be instances where the set is convex, even though the
answer of this function is `false`.
"""
function isconvextype(::Type{<:LazySet}) end
