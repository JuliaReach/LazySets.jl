import IntervalArithmetic: ⊂

"""
    ⊂(X::LazySet{N}, Y::LazySet, [witness]::Bool=false) where {N}

Strict inclusion check.

### Input

- `X`       -- first set
- `Y`       -- second set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true` iff ``X ⊂ Y``
* If `witness` option is activated:
  * `(true, v)` iff ``X ⊂ Y`` and ``v ∈ Y \\setminus X``
  * `(false, [])` iff not ``X ⊂ Y``

### Algorithm

We check inclusion of `X` in `Y` and then check inclusion of `Y` in `X`:

```math
X ⊂ Y \\Leftrightarrow X ⊆ Y \\land ¬ (Y ⊆ X)
```
"""
function ⊂(X::LazySet{N}, Y::LazySet, witness::Bool=false) where {N}
    if witness
        res, w = ⊆(X, Y, witness)
        if res
            res, w = ⊆(Y, X, witness)
            if res
                return (false, N[])
            else
                return (true, w)
            end
        else
            return (res, N[])
        end
    end
    return (X ⊆ Y) && !(Y ⊆ X)
end
