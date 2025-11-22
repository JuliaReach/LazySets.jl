"""
# Extended help

    linear_combination(P1::SimpleSparsePolynomialZonotope,
                       P2::SimpleSparsePolynomialZonotope)

### Notes

This method implements the algorithm described in [Kochdumper21a; Proposition 3.1.25](@citet).
"""
@validate function linear_combination(P1::SimpleSparsePolynomialZonotope,
                                      P2::SimpleSparsePolynomialZonotope)
    c1, c2 = center(P1), center(P2)
    G1, G2 = genmat(P1), genmat(P2)
    E1, E2 = expmat(P1), expmat(P2)

    return _linear_combination(c1, c2, G1, G2, E1, E2)
end

function _linear_combination(c1, c2, G1, G2, E1, E2)
    h1, h2 = size(G1, 2), size(G2, 2)
    p1, p2 = size(E1, 1), size(E2, 1)

    c = (c1 + c2) / 2
    G = hcat(c1 - c2, G1, G1, G2, -G2) / 2
    E = hcat(vcat(zeros(Int, p1 + p2), 1),
             vcat(E1, zeros(Int, p2 + 1, h1)),
             vcat(E1, zeros(Int, p2, h1), ones(Int, 1, h1)),
             vcat(zeros(Int, p1, h2), E2, zeros(Int, 1, h2)),
             vcat(zeros(Int, p1, h2), E2, ones(Int, 1, h2)))

    return SimpleSparsePolynomialZonotope(c, G, E)
end
