"""
    minkowski_sum(P1::SparsePolynomialZonotope, P2::SparsePolynomialZonotope)

Compute the Minkowski sum of two sparse polyomial zonotopes.

### Input

- `P1` -- sparse polynomial zonotope
- `P2` -- sparse polynomial zonotope

### Output

The Minkowski sum of `P1` and `P2`.

### Algorithm

See Proposition 3.1.19 in [1].

[1] Kochdumper. *Extensions of polynomial zonotopes and their application to
    verification of cyber-physical systems.* PhD diss., TU Munich, 2022.
"""
function minkowski_sum(P1::SparsePolynomialZonotope,
                       P2::SparsePolynomialZonotope)
    c = center(P1) + center(P2)
    G = hcat(genmat_dep(P1), genmat_dep(P2))
    GI = hcat(genmat_indep(P1), genmat_indep(P2))
    E = cat(expmat(P1), expmat(P2); dims=(1, 2))
    return SparsePolynomialZonotope(c, G, GI, E)
end
