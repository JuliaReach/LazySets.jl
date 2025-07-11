"""
    indexvector(MZ::MatrixZonotope)

Return the  index vector of the matrix zonotope `MZ`.

### Input

- `MZ` -- matrix zonotope set

### Output

A vector of unique positive integers representing each generator
"""
function indexvector(MZ::MatrixZonotope)
    return MZ.idx
end