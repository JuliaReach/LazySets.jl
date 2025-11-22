"""
# Extended help

    convex_hull(P::SimpleSparsePolynomialZonotope)

### Output

The tightest convex simple sparse polynomial zonotope containing `P`.
"""
function convex_hull(P::SimpleSparsePolynomialZonotope)
    return linear_combination(P, P)
end

# [Kochdumper21a; Proposition 3.1.27](@citet)
@validate function convex_hull(P1::SimpleSparsePolynomialZonotope,
                               P2::SimpleSparsePolynomialZonotope)
    c1, c2 = center(P1), center(P2)
    G1, G2 = genmat(P1), genmat(P2)
    G1bar = hcat(G1, G1, G1, -G1) / 2
    G2bar = hcat(G2, G2, G2, -G2) / 2
    E1, E2 = expmat(P1), expmat(P2)
    E1bar = _Ebar(E1)
    E2bar = _Ebar(E2)

    return _linear_combination(c1, c2, G1bar, G2bar, E1bar, E2bar)
end

function _Ebar(E)
    p, h = size(E)
    return hcat(vcat(E, zeros(Int, p + 1, h)),
                vcat(E, zeros(Int, p, h), ones(Int, 1, h)),
                vcat(zeros(Int, p, h), E, zeros(Int, 1, h)),
                vcat(zeros(Int, p, h), E, ones(Int, 1, h)))
end
