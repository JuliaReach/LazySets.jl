"""
    isdisjoint(X::LazySet, Y::LazySet, [witness]::Bool=false)

Check whether two sets are disjoint (i.e., do not intersect), and optionally
compute a witness.

### Input

- `X`       -- set
- `Y`       -- set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If the `witness` option is deactivated: `true` iff ``X ∩ Y = ∅``
* If the `witness` option is activated:
  * `(true, [])` iff ``X ∩ Y = ∅``
  * `(false, v)` iff ``X ∩ Y ≠ ∅`` for some ``v ∈ X ∩ Y``

### Notes

The convenience alias `is_intersection_empty` is also available.
"""
function isdisjoint(::LazySet, ::LazySet, ::Bool=false) end

const is_intersection_empty = isdisjoint
