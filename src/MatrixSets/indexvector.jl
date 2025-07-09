"""
    indexvector(MZ::MatrixZonotope)

Return the `indexVector` of the matrix zonotope `Z`.

### Input

- `MZ` -- matrix zonotope set

### Output

A vector of unique positive integers representing each generator
"""
function indexvector(MZ::MatrixZonotope)
    return MZ.idx
end