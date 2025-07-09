"""
    center(MZ::MatrixZonotope)

Return the center matrix of the matrix zonotope `Z`.

### Input

- `MZ` -- matrix zonotope set

### Output

The centre matrix of `MZ`.
"""
center(MZ::MatrixZonotope) = MZ.A0 
