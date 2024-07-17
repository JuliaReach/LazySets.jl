"""
    isbounded(hs::HalfSpace)

Check whether a half-space is bounded.

### Input

- `hs` -- half-space

### Output

`false`.
"""
function isbounded(::HalfSpace)
    return false
end
