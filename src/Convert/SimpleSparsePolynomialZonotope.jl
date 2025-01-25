"""
    convert(::Type{SimpleSparsePolynomialZonotope}, Z::AbstractZonotope)

Convert a zonotope to a simple sparse polynomial zonotope.

### Input

- `SimpleSparsePolynomialZonotope` -- target type
- `Z`                              -- zonotopic set

### Output

A simple sparse polynomial zonotope.

### Algorithm

This method implements [KochdumperA21; Proposition 3](@citet).
"""
function Base.convert(::Type{SimpleSparsePolynomialZonotope}, Z::AbstractZonotope)
    c = center(Z)
    G = genmat(Z)
    n = ngens(Z)
    E = Matrix(1 * I, n, n)
    return SimpleSparsePolynomialZonotope(c, G, E)
end

"""
    convert(::Type{SimpleSparsePolynomialZonotope}, SPZ::SparsePolynomialZonotope)

Convert a sparse polynomial zonotope to simple sparse polynomial zonotope.

### Input

- `SimpleSparsePolynomialZonotope` -- target type
- `SPZ`                            -- sparse polynomial zonotope

### Output

A simple sparse polynomial zonotope.

### Algorithm

The method implements [Kochdumper21a; Proposition 3.1.4](@citet).
"""
function Base.convert(::Type{SimpleSparsePolynomialZonotope}, SPZ::SparsePolynomialZonotope)
    c = center(SPZ)
    G = hcat(genmat_dep(SPZ), genmat_indep(SPZ))
    n = ngens_indep(SPZ)
    E = cat(expmat(SPZ), Matrix(1 * I, n, n); dims=(1, 2))
    return SimpleSparsePolynomialZonotope(c, G, E)
end
