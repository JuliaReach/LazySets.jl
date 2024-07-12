"""
    translate(S::SparsePolynomialZonotope, v::AbstractVector)

Translate (i.e., shift) a sparse polynomial zonotope by a given vector.

### Input

- `S` -- sparse polynomial zonotope
- `v` -- translation vector

### Output

A translated sparse polynomial zonotope.
"""
function translate(S::SparsePolynomialZonotope, v::AbstractVector)
    c = center(S) + v
    return SparsePolynomialZonotope(c, genmat_dep(S), genmat_indep(S), expmat(S))
end
