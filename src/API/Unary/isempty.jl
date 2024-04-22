"""
    isempty(X::LazySet, witness::Bool=false)

Check whether a set is empty.

### Input

- `X`       -- set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If the `witness` option is deactivated: `true` iff ``X = ∅``
* If the `witness` option is activated:
  * `(true, [])` iff ``X = ∅``
  * `(false, v)` iff ``X ≠ ∅`` for some ``v ∈ X``
"""
function isempty(::LazySet, ::Bool=false) end
