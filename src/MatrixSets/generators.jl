"""
    generators(MZ::MatrixZonotope)

Return the generators of a matrix zonotope.

### Input

- `MZ` -- matrix zonotope set

### Output

The generators of `MZ`.
"""
generators(MZ::MatrixZonotope) = MZ.Ai 
