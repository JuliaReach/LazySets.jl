"""
    linear_map(M::AbstractMatrix, P::DensePolynomialZonotope)

Return the linear map of a polynomial zonotope.

### Input

- `M` -- matrix
- `P` -- polynomial zonotope

## Output

A polynomial zonotope.

### Algorithm

The result's starting point and generators are those of `P` multiplied by the
matrix `M`.
"""
function linear_map(M::AbstractMatrix, P::DensePolynomialZonotope)
    c = M * P.c
    E = [M * Ei for Ei in P.E]
    F = [M * Fi for Fi in P.F]
    G = M * P.G
    return DensePolynomialZonotope(c, E, F, G)
end
