"""
    isuniversal(L::Line2D, [witness]::Bool=false)

Check whether a 2D line is universal.

### Input

- `L`       -- 2D line
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `false`
* If `witness` option is activated: `(false, v)` where ``v ∉ L``

### Algorithm

Witness production falls back to `isuniversal(::Hyperplane)`.
"""
function isuniversal(L::Line2D, witness::Bool=false)
    if witness
        v = _non_element_halfspace(L.a, L.b)
        return (false, v)
    else
        return false
    end
end
