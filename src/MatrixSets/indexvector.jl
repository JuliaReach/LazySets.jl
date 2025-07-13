"""
    indexvector(MZ::MatrixZonotope)

Return the index vector of a matrix zonotope.

### Input

- `MZ` -- matrix zonotope

### Output

A vector of unique positive integers representing each generator.
"""
function indexvector(MZ::MatrixZonotope)
    return MZ.idx
end
