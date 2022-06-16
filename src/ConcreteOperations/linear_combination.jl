export linear_combination

"""
    linear_combination(P1::SimpleSparsePolynomialZonotope, P2::SimpleSparsePolynomialZonotope)

compute the linear combination of simple sparse polynomial zonotopes `P1` and `P2`.

### Input

- `P1` -- simple sparse polynomial zonotope
- `P2` -- simple sparse polynomial zonotope

### Output

linear combination of `P1` and `P2`.

### Notes

The linear combination of two sets ``P₁`` and ``P₂`` is defined as

```math
\\{λp₁ + (1-λ)p₂ | p₁ ∈ P₁, p₂ ∈ P₂, λ ∈ [0, 1]\\}
```
"""
function linear_combination(P1::SimpleSparsePolynomialZonotope, P2::SimpleSparsePolynomialZonotope)
    c = 0.5 * (center(P1) + center(P2))
    G = 0.5 * hcat(center(P1) - center(P2), genmat(P1), genmat(P1), genmat(P2), -genmat(P2))
    E1, E2 = expmat(P1), expmat(P2)
    et = promote_type(eltype(E1), eltype(E2))
    E = hcat(vcat(zeros(et, nparams(P1) + nparams(P2)), one(et)),
             vcat(E1, zeros(et, nparams(P2) + 1, ngens(P1))),
             vcat(E2, zeros(et, nparams(P2), ngens(P1)), ones(et, 1, ngens(P1))),
             vcat(zeros(et, nparams(P1), ngens(P2)), E2, zeros(et, 1, ngens(P2))),
             vcat(zeros(et, nparams(P1), ngens(P2)), E2, ones(et, 1, ngens(P2))))

    return SimpleSparsePolynomialZonotope(c, G, E)
end
