"""
    order(MZ::MatrixZonotope)

Return the order of a matrix zonotope.

### Input

- `MZ` -- matrix zonotope

### Output

A rational number representing the order of the matrix zonotope.

### Notes

The order of a matrix zonotope is defined as the quotient of its number of
generators and the product of its dimensions. Alternatively it can be thought
as the order of the zonotopic set constructed from the vectorization of the center 
and the generators.
"""
function order(MZ::MatrixZonotope)
    m, n = size(MZ)
    return ngens(MZ) // (m * n)
end