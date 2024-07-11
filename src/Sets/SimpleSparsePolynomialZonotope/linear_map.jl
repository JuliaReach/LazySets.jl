"""
    linear_map(M::AbstractMatrix, P::SimpleSparsePolynomialZonotope)

Apply the linear map `M` to a simple sparse polynomial zonotope.

### Input

- `M` -- matrix
- `P` -- simple sparse polynomial zonotope

### Output

The set resulting from applying the linear map `M` to `P`.
"""
function linear_map(M::AbstractMatrix, P::SSPZ)
    return SimpleSparsePolynomialZonotope(M * center(P), M * genmat(P), expmat(P))
end
