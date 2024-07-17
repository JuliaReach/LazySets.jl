"""
    constrained_dimensions(hs::HalfSpace)

Return the indices in which a half-space is constrained.

### Input

- `hs` -- half-space

### Output

A vector of ascending indices `i` such that the half-space is constrained in
dimension `i`.

### Examples

A 2D half-space with constraint ``x_1 ≥ 0`` is only constrained in dimension 1.
"""
function constrained_dimensions(hs::HalfSpace)
    return nonzero_indices(hs.a)
end
