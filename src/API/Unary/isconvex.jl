"""
    isconvex(X::LazySet)

Check whether a set is convex.

### Input

- `X` -- set

### Output

`true` iff the set is convex.

### Notes

The answer is conservative, i.e., may sometimes be `false` even if the set is
convex.

See also [`isconvextype(::Type{<:LazySet})`](@ref).
"""
function isconvex(::LazySet) end
