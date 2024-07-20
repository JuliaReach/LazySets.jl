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

"""
    convex_hull(P1::SimpleSparsePolynomialZonotope,
                P2::SimpleSparsePolynomialZonotope)

Compute the convex hull of two simple sparse polynomial zonotopes.

### Input

- `P1` : simple sparse polynomial zonotopes
- `P2` : simple sparse polynomial zonotopes

### Output

Tightest convex simple sparse polynomial zonotope containing `P1` and `P2`.
"""
function convex_hull(P1::SimpleSparsePolynomialZonotope, P2::SimpleSparsePolynomialZonotope)
    return linear_combination(linear_combination(P1, P1), linear_combination(P2, P2))
end
