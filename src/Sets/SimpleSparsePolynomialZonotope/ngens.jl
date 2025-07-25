"""
# Extended help

    ngens(P::SimpleSparsePolynomialZonotope)

### Notes

This number corresponds to the number of monomials in the polynomial
representation of `P`.
"""
ngens(P::SimpleSparsePolynomialZonotope) = ngens_dep(P)
