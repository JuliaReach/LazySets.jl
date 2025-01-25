"""
# Extended help

    cartesian_product(P1::SimpleSparsePolynomialZonotope,
                      P2::SimpleSparsePolynomialZonotope)

### Algorithm

This method implements [Kochdumper21a; Proposition 3.1.22](@citet).
"""
function cartesian_product(P1::SimpleSparsePolynomialZonotope,
                           P2::SimpleSparsePolynomialZonotope)
    c = vcat(center(P1), center(P2))
    G = cat(genmat(P1), genmat(P2); dims=(1, 2))
    E = cat(expmat(P1), expmat(P2); dims=(1, 2))
    return SimpleSparsePolynomialZonotope(c, G, E)
end
