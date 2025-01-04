"""
# Extended help

    convex_hull(P::SimpleSparsePolynomialZonotope)

### Output

The tightest convex simple sparse polynomial zonotope containing `P`.
"""
function convex_hull(P::SimpleSparsePolynomialZonotope)
    return linear_combination(P, P)
end
