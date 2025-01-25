"""
# Extended help

    exact_sum(P1::SparsePolynomialZonotope, P2::SparsePolynomialZonotope)

### Algorithm

This method implements [Kochdumper21a; Proposition 3.1.20](@citet).
"""
function exact_sum(P1::SparsePolynomialZonotope, P2::SparsePolynomialZonotope)
    if indexvector(P1) != indexvector(P2)
        throw(ArgumentError("the exact sum is currently only implemented for " *
                            "sparse polynomial zonotopes with the same index vector"))
    end

    c = center(P1) + center(P2)
    G = hcat(genmat_dep(P1), genmat_dep(P2))
    GI = hcat(genmat_indep(P1), genmat_indep(P2))
    E = hcat(expmat(P1), expmat(P2))
    idx = indexvector(P1)

    return SparsePolynomialZonotope(c, G, GI, E, idx)
end
