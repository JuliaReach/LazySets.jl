"""
    ⊆(x::Interval, y::Interval, [witness]::Bool=false)

Check whether an interval is contained in another interval.

### Input

- `x`       -- inner interval
- `y`       -- outer interval
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

`true` iff ``x ⊆ y``.
"""
function ⊆(x::Interval, y::Interval, witness::Bool=false)
    if min(y) > min(x)
        return witness ? (false, low(x)) : false
    elseif max(x) > max(y)
        return witness ? (false, high(x)) : false
    end
    return _witness_result_empty(witness, true, x, y)
end
