"""
    isequivalent(X::LazySet, Y::LazySet)

Check whether two sets are equivalent, i.e., represent the same set of points.

### Input

- `X` -- set
- `Y` -- set

### Output

`true` iff `X` is equivalent to `Y` (up to numerical precision).
"""
function isequivalent(::LazySet, ::LazySet) end
