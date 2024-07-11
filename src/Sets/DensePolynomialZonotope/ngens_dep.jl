"""
    ngens_dep(P::DensePolynomialZonotope)

Return the number of dependent generators of a polynomial zonotope.

### Input

- `P` -- polynomial zonotope

### Output

The number of dependent generators of `P`.
"""
function ngens_dep(P::DensePolynomialZonotope)
    η = polynomial_order(P)  # polynomial order
    p = size(P.E[1], 2)  # number of dependent factors
    return sum(i -> binomial(p + i - 1, i), 1:η)
end
