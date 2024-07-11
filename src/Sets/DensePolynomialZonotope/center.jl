"""
    center(P::DensePolynomialZonotope)

Return the center of a polynomial zonotope.

### Input

- `P` -- polynomial zonotope

### Output

The center of `P`.
"""
center(P::DensePolynomialZonotope) = P.c
