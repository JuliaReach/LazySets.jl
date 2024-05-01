"""
    isuniversal(X::LazySet, witness::Bool=false)

Check whether a set is universal.

### Input

- `X`       -- set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If the `witness` option is deactivated: `true` iff ``X = ℝ^n``
* If the `witness` option is activated:
  * `(true, [])` iff ``X = ℝ^n``
  * `(false, v)` iff ``X ≠ ℝ^n`` for some ``v ∉ X``
"""
function isuniversal(::LazySet, ::Bool=false) end
