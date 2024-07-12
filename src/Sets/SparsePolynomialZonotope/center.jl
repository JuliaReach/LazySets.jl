"""
    center(P::SparsePolynomialZonotope)

Return the center of a sparse polynomial zonotope.

### Input

- `P` -- sparse polynomial zonotope

### Output

The center.
"""
center(P::SPZ) = P.c
