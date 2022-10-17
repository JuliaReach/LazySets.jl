export exact_sum, ⊞

"""
    exact_sum(P1::SparsePolynomialZonotope, P2::SparsePolynomialZonotope)

Compute the exact sum of sparse polyomial zonotopes ``P₁`` and ``P₂``.

### Input

- `P1` -- sparse polynomial zonotope
- `P2` -- sparse polynomial zonotope

### Output

A `SparsePolynomialZonotope` representing the exact sum ``P₁ ⊞ P₂``.
"""
function exact_sum(P1::SparsePolynomialZonotope, P2::SparsePolynomialZonotope)

    indexvector(P1) == indexvector(P2) || throw(ArgumentError("the exact sum " *
        "is currently only implemented for sparse polynomial zonotopes with " *
        "the same index vector"))

    c = center(P1) + center(P2)
    G = hcat(genmat_dep(P1), genmat_dep(P2))
    GI = hcat(genmat_indep(P1), genmat_indep(P2))
    E = hcat(expmat(P1), expmat(P2))
    idx = indexvector(P1)

    return SparsePolynomialZonotope(c, G, GI, E, idx)
end

"""
    ⊞(X::LazySet, Y::LazySet)

Unicode alias constructor for the (concrete) `exact_sum` function.

### Notes

Write `\\boxplus[TAB]` to enter this symbol.
"""
const ⊞ = exact_sum
