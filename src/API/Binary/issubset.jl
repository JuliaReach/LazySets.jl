"""
    ⊆(X::LazySet, Y::LazySet, [witness]::Bool=false)

Check whether a set is a subset of another set, and optionally compute a witness.

### Input

- `X`       -- set
- `Y`       -- set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If the `witness` option is deactivated: `true` iff ``X ⊆ Y``
* If the `witness` option is activated:
  * `(true, [])` iff ``X ⊆ Y``
  * `(false, v)` iff ``X ⊈ Y`` for some ``v ∈ X ∖ Y``

### Notes

The convenience alias `issubset` is also available.
"""
function ⊆(::LazySet, ::LazySet, ::Bool=false) end
