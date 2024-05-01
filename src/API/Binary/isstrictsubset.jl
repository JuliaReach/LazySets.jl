"""
    ⊂(X::LazySet, Y::LazySet, [witness]::Bool=false)

Check whether a set is a strict subset of another set, and optionally compute a
witness.

### Input

- `X`       -- set
- `Y`       -- set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If the `witness` option is deactivated: `true` iff ``X ⊂ Y``
* If the `witness` option is activated:
  * `(true, v)` iff ``X ⊂ Y`` for some ``v ∈ Y ∖ X``
  * `(false, [])` iff ``X ⊂ Y`` does not hold

### Notes

Strict inclusion is sometimes written as `⊊`. The following identity holds:

```math
    X ⊂ Y ⇔ X ⊆ Y ∧ Y ⊈ X
```

### Algorithm

The default implementation first checks inclusion of `X` in `Y` and then checks
noninclusion of `Y` in `X`:
"""
function ⊂(::LazySet, ::LazySet, ::Bool=false) end
