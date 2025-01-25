"""
# Extended help

    minkowski_sum(P1::SparsePolynomialZonotope, P2::SparsePolynomialZonotope)

### Algorithm

See [Kochdumper21a; Proposition 3.1.19](@citet).
"""
function minkowski_sum(P1::SparsePolynomialZonotope,
                       P2::SparsePolynomialZonotope)
    c = center(P1) + center(P2)
    G = hcat(genmat_dep(P1), genmat_dep(P2))
    GI = hcat(genmat_indep(P1), genmat_indep(P2))
    E = cat(expmat(P1), expmat(P2); dims=(1, 2))
    return SparsePolynomialZonotope(c, G, GI, E)
end
