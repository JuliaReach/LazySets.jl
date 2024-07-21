"""
    minkowski_sum(P1::SimpleSparsePolynomialZonotope,
                  P2::SimpleSparsePolynomialZonotope)

Compute the Minkowski sum of two simple sparse polynomial zonotopes.

### Input

- `P1` -- simple sparse polynomial zonotope
- `P2` -- simple sparse polynomial zonotope

### Output

The Minkowski sum of `P1` and `P2`.
"""
function minkowski_sum(P1::SimpleSparsePolynomialZonotope,
                       P2::SimpleSparsePolynomialZonotope)
    c = center(P1) + center(P2)
    G = hcat(genmat(P1), genmat(P2))
    E = cat(expmat(P1), expmat(P2); dims=(1, 2))
    return SimpleSparsePolynomialZonotope(c, G, E)
end
