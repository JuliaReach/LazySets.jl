"""
    convex_hull(P::SimpleSparsePolynomialZonotope)

Compute the convex hull of a simple sparse polynomial zonotope.

### Input

- `P` -- simple sparse polynomial zonotope

### Output

The tightest convex simple sparse polynomial zonotope containing `P`.
"""
function convex_hull(P::SimpleSparsePolynomialZonotope)
    return linear_combination(P, P)
end
