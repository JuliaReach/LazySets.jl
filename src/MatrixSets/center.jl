"""
    center(MZ::MatrixZonotope)

Return the center matrix of a matrix zonotope.

### Input

- `MZ` -- matrix zonotope

### Output

The center matrix of `MZ`.
"""
function center(MZ::MatrixZonotope)
    return MZ.A0
end
