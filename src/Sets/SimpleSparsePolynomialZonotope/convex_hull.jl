"""
# Extended help

    convex_hull(P::SimpleSparsePolynomialZonotope)

### Output

The tightest convex simple sparse polynomial zonotope containing `P`.
"""
function convex_hull(P::SimpleSparsePolynomialZonotope)
    return linear_combination(P, P)
end

"""
# Extended help

    convex_hull(P1::SimpleSparsePolynomialZonotope,
                P2::SimpleSparsePolynomialZonotope)

### Output

The tightest convex simple sparse polynomial zonotope containing `P1` and `P2`.
"""
@validate function convex_hull(P1::SimpleSparsePolynomialZonotope,
                               P2::SimpleSparsePolynomialZonotope)
    return linear_combination(linear_combination(P1, P1), linear_combination(P2, P2))
end
