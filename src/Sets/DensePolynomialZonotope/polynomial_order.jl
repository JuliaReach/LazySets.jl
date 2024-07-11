"""
    polynomial_order(P::DensePolynomialZonotope)

Polynomial order of a polynomial zonotope.

### Input

- `P` -- polynomial zonotope

## Output

The polynomial order, defined as the maximal power of the scale factors ``β_i``.
It is usually denoted ``η``.
"""
polynomial_order(P::DensePolynomialZonotope) = length(P.E)
