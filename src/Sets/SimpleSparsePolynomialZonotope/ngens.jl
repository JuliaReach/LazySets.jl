"""
    ngens(P::SimpleSparsePolynomialZonotope)

Return the number of generators of a simple sparse polynomial zonotope.

### Input

- `P` -- simple sparse polynomial zonotope

### Output

The number of generators of `P`.

### Notes

This number corresponds to the number of monomials in the polynomial
representation of `P`.
"""
ngens(P::SSPZ) = ngens_dep(P)
