"""
    convert(::Type{SimpleSparsePolynomialZonotope}, SPZ::AbstractSparsePolynomialZonotope)

Convert a sparse polynomial zonotope to simple sparse polynomial zonotope.

### Input

- `SimpleSparsePolynomialZonotope` -- target type
- `SPZ`                            -- sparse polynomial zonotope

### Output

A simple sparse polynomial zonotope.

### Algorithm

The method implements [Kochdumper21a; Proposition 3.1.4](@citet).
"""
function convert(::Type{SimpleSparsePolynomialZonotope}, SPZ::AbstractSparsePolynomialZonotope)
    c = center(SPZ)
    G = hcat(genmat_dep(SPZ), genmat_indep(SPZ))
    n = ngens_indep(SPZ)
    E = cat(expmat(SPZ), Matrix(1 * I, n, n); dims=(1, 2))
    return SimpleSparsePolynomialZonotope(c, G, E)
end
