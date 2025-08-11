"""
    remove_redundant_generators(MZ::MatrixZonotope)

Remove redundant generators from a matrix zonotope.

# Input

- `MZ` -- a matrix zonotope

# Output

A new matrix zonotope with fewer generators, or the same matrix zonotope 
if no generator could be removed.

# Algorithm 

This function first vectorizes the matrix zonotope into a standard zonotope, 
removes redundant generators from the resulting zonotope, and then converts it back 
to a matrix zonotope of the original dimensions.

# Extended help

    remove_redundant_generators(Z::Zonotope)
"""
function remove_redundant_generators(MZ::MatrixZonotope)
    Z = vectorize(MZ)
    Zred = remove_redundant_generators(Z)

    dim = size(MZ)
    return matrixize(Zred, dim)
end
