"""
    expmat(P::SparsePolynomialZonotope)

Return the matrix of exponents of the sparse polynomial zonotope.

### Input

- `P` -- sparse polynomial zonotope

### Output

The matrix of exponents, where each column is a multidegree.

### Notes

In the exponent matrix, each row corresponds to a parameter (``αₖ`` in the
definition) and each column to a monomial.
"""
expmat(P::SPZ) = P.E
