"""
    isuniversal(U::Universe{N}, [witness]::Bool=false) where {N}

Check whether a universe is universal.

### Input

- `U`       -- universe
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `true`
* If `witness` option is activated: `(true, [])`
"""
function isuniversal(U::Universe{N}, witness::Bool=false) where {N}
    return witness ? (true, N[]) : true
end
