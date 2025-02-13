"""
    convex_hull(X::LazySet, Y::LazySet)

Compute the convex hull of (the union of) two sets.

### Input

- `X` -- set
- `Y` -- set

### Output

A set representing the convex hull of ``X ∪ Y``.

### Notes

See [`convex_hull(::LazySet)`](@ref) for the convex hull of a single set.
"""
function convex_hull(::LazySet, ::LazySet) end
