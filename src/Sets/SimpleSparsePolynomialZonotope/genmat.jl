"""
    genmat(P::SimpleSparsePolynomialZonotope)

Return the matrix of generators of a simple sparse polynomial zonotope.

### Input

- `P` -- simple sparse polynomial zonotope

### Output

The matrix of generators of `P`.
"""
genmat(P::SSPZ) = P.G
