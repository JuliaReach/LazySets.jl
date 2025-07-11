"""
    center(MZ::MatrixZonotope)

Return the center matrix of the matrix zonotope `Z`.

### Input

- `MZ` -- matrix zonotope set

### Output

The center matrix of `MZ`.
"""
center(MZ::MatrixZonotope) = MZ.A0 
