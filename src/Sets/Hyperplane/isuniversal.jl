"""
    isuniversal(H::Hyperplane, [witness]::Bool=false)

Check whether a hyperplane is universal.

### Input

- `P`       -- hyperplane
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `false`
* If `witness` option is activated: `(false, v)` where ``v ∉ P``

### Algorithm

A witness is produced by adding the normal vector to an element on the
hyperplane.
"""
function isuniversal(H::Hyperplane, witness::Bool=false)
    if witness
        v = _non_element_halfspace(H.a, H.b)
        return (false, v)
    else
        return false
    end
end
