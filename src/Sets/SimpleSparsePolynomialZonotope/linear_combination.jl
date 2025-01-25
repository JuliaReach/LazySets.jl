"""
# Extended help

    linear_combination(P1::SimpleSparsePolynomialZonotope,
                       P2::SimpleSparsePolynomialZonotope)

### Notes

This method implements the algorithm described in [Kochdumper21a; Proposition 3.1.25](@citet).
"""
function linear_combination(P1::SimpleSparsePolynomialZonotope,
                            P2::SimpleSparsePolynomialZonotope)
    c1, c2 = center(P1), center(P2)
    G1, G2 = genmat(P1), genmat(P2)
    E1, E2 = expmat(P1), expmat(P2)

    c = (c1 + c2) / 2
    G = hcat(c1 - c2, G1, G1, G2, -G2) / 2

    N = promote_type(eltype(E1), eltype(E2))
    E = hcat(vcat(zeros(N, nparams(P1) + nparams(P2)), one(N)),
             vcat(E1, zeros(N, nparams(P2) + 1, ngens(P1))),
             vcat(E1, zeros(N, nparams(P2), ngens(P1)), ones(N, 1, ngens(P1))),
             vcat(zeros(N, nparams(P1), ngens(P2)), E2, zeros(N, 1, ngens(P2))),
             vcat(zeros(N, nparams(P1), ngens(P2)), E2, ones(N, 1, ngens(P2))))

    return SimpleSparsePolynomialZonotope(c, G, E)
end
