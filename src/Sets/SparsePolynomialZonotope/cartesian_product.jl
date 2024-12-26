"""
# Extended help

    cartesian_product(P1::SparsePolynomialZonotope, P2::SparsePolynomialZonotope)

### Algorithm

This method implements Proposition 3.1.22 in [1].

[1] Kochdumper, Niklas. *Extensions of polynomial zonotopes and their application
    to verification of cyber-physical systems.* PhD diss., Technische Universität
    München, 2022.
"""
function cartesian_product(P1::SparsePolynomialZonotope, P2::SparsePolynomialZonotope)
    c = vcat(center(P1), center(P2))
    G = cat(genmat_dep(P1), genmat_dep(P2); dims=(1, 2))
    GI = cat(genmat_indep(P1), genmat_indep(P2); dims=(1, 2))
    E = cat(expmat(P1), expmat(P2); dims=(1, 2))
    return SparsePolynomialZonotope(c, G, GI, E)
end
