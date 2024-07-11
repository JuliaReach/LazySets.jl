"""
    ngens_indep(P::DensePolynomialZonotope)

Return the number of independent generators of a polynomial zonotope.

### Input

- `P` -- polynomial zonotope

### Output

The number of independent generators of `P`.
"""
ngens_indep(P::DensePolynomialZonotope) = size(P.G, 2)
