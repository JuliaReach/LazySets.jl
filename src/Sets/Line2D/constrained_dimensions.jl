"""
    constrained_dimensions(L::Line2D)

Return the indices in which a 2D line is constrained.

### Input

- `L` -- 2D line

### Output

A vector of ascending indices `i` such that the line is constrained in dimension
`i`.

### Examples

A line with constraint ``x_i = 0`` (``i ∈ \\{1, 2\\}``) is only constrained in
dimension ``i``.
"""
function constrained_dimensions(L::Line2D)
    return nonzero_indices(L.a)
end
