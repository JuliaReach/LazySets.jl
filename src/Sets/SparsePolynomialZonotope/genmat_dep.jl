"""
    genmat_dep(P::SparsePolynomialZonotope)

Return the matrix of dependent generators of a sparse polynomial zonotope.

### Input

- `P` -- sparse polynomial zonotope

### Output

The matrix of dependent generators.
"""
genmat_dep(P::SPZ) = P.G
