"""
    dim(hs::HalfSpace)

Return the dimension of a half-space.

### Input

- `hs` -- half-space

### Output

The ambient dimension of the half-space.
"""
function dim(hs::HalfSpace)
    return length(hs.a)
end
