"""
    isbounded(B::Bloating)

Determine whether a bloated set is bounded.

### Input

- `B` -- bloated set

### Output

`true` iff the wrapped set is bounded.
"""
function isbounded(B::Bloating)
    return isbounded(B.X)
end
