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
function convert(::Type{SimpleSparsePolynomialZonotope}, Z::AbstractZonotope)
    c = center(Z)
    G = genmat(Z)
    n = ngens(Z)
    E = Matrix(1 * I, n, n)
    return SimpleSparsePolynomialZonotope(c, G, E)
end
