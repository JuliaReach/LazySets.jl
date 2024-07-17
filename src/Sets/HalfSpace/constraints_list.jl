"""
    constraints_list(hs::HalfSpace)

Return the list of constraints of a half-space.

### Input

- `hs` -- half-space

### Output

A singleton list containing the half-space.
"""
function constraints_list(hs::HalfSpace)
    return [hs]
end
