"""
	generators(MZ::MatrixZonotope)

Return an iterator over the generators of a matrix zonotope.

### Input

- `MZ` -- matrix zonotope set

### Output

An iterator over the generators of `MZ`.
"""
generators(MZ::MatrixZonotope) = MZ.Ai
