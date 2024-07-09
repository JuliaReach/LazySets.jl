"""
    shape_matrix(E::Ellipsoid)

Return the shape matrix of the ellipsoid.

### Input

- `E` -- ellipsoid

### Output

The shape matrix of the ellipsoid.
"""
function shape_matrix(E::Ellipsoid)
    return E.shape_matrix
end
