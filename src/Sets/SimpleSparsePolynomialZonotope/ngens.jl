"""
# Extended help

    ngens(P::SimpleSparsePolynomialZonotope)

### Notes

This number corresponds to the number of monomials in the polynomial
representation of `P`.
"""
function ngens(P::SimpleSparsePolynomialZonotope)
    return ngens_dep(P)
end
