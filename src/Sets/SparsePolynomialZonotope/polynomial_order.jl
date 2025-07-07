"""
    polynomial_order(P::SparsePolynomialZonotope)

Return the polynomial order of a sparse polynomial zonotope.

### Input

- `P` -- sparse polynomial zonotope

### Output

The polynomial order.

### Notes

The polynomial order is the maximum sum of all monomials' parameter exponents.
"""
function polynomial_order(P::SparsePolynomialZonotope)
    return maximum(sum, eachcol(expmat(P)))
end
