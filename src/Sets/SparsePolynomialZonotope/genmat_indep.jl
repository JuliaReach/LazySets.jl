"""
    genmat_indep(P::SparsePolynomialZonotope)

Return the matrix of independent generators of a sparse polynomial zonotope.

### Input

- `P` -- sparse polynomial zonotope

### Output

The matrix of independent generators.
"""
genmat_indep(P::SparsePolynomialZonotope) = P.GI
