"""
    center(P::SimpleSparsePolynomialZonotope)

Return the center of a simple sparse polynomial zonotope.

### Input

- `P` -- simple sparse polynomial zonotope

### Output

The center of `P`.
"""
center(P::SSPZ) = P.c
