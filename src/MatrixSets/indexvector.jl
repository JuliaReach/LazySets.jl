"""
	indexvector(MZ::MatrixZonotope)

Return the `indexVector` of the matrix zonotope `Z`.

### Input

- `MZ` -- matrix zonotope set

### Output

The index vector contains positive integers for the generators.
"""
function indexvector(MZ::MatrixZonotope)
    return MZ.idx
end