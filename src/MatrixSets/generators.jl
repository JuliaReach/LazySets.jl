"""
    generators(MZ::MatrixZonotope)

Return the generators of a matrix zonotope.

### Input

- `MZ` -- matrix zonotope

### Output

The generators of `MZ`.
"""
function generators(MZ::MatrixZonotope)
    return MZ.Ai
end
