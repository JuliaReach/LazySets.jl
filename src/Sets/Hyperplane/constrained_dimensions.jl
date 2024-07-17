"""
    constrained_dimensions(H::Hyperplane)

Return the dimensions in which a hyperplane is constrained.

### Input

- `H` -- hyperplane

### Output

A vector of ascending indices `i` such that the hyperplane is constrained in
dimension `i`.

### Examples

A 2D hyperplane with constraint ``x_1 = 0`` is constrained in dimension 1 only.
"""
function constrained_dimensions(H::Hyperplane)
    return nonzero_indices(H.a)
end
