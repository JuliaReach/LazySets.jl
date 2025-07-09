"""
    ngens(MZ::MatrixZonotope)

Return the number of generators of a matrix zonotope.

### Input

- `MZ` -- matrix zonotope 

### Output

An integer representing the number of generators.
"""
function ngens(MZ::MatrixZonotope)
    return length(generators(MZ))
end
