"""
    isuniversal(∅::EmptySet{N}, [witness]::Bool=false) where {N}

Check whether an empty set is universal.

### Input

- `∅`       -- empty set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `false`
* If `witness` option is activated: `(false, v)` where ``v ∉ S``
"""
function isuniversal(∅::EmptySet{N}, witness::Bool=false) where {N}
    if witness
        return (false, zeros(N, dim(∅)))
    else
        return false
    end
end
