"""
    cartesian_product(P1::SimpleSparsePolynomialZonotope,
                      P2::SimpleSparsePolynomialZonotope)

Compute the Cartesian product of two simple sparse polynomial zonotopes.

### Input

- `P1` -- simple sparse polynomial zonotope
- `P2` -- simple sparse polynomial zonotope

### Output

The Cartesian product of `P1` and `P2`.

### Algorithm

This method implements Proposition 3.1.22 in [1].

[1] Kochdumper, Niklas. *Extensions of polynomial zonotopes and their application
    to verification of cyber-physical systems.* PhD diss., Technische Universität
    München, 2022.
"""
function cartesian_product(P1::SimpleSparsePolynomialZonotope,
                           P2::SimpleSparsePolynomialZonotope)
    c = vcat(center(P1), center(P2))
    G = cat(genmat(P1), genmat(P2); dims=(1, 2))
    E = cat(expmat(P1), expmat(P2); dims=(1, 2))
    return SimpleSparsePolynomialZonotope(c, G, E)
end
