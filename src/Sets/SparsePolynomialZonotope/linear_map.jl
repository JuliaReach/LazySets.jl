"""
    linear_map(M::AbstractMatrix, P::SparsePolynomialZonotope)

Apply a linear map to a sparse polynomial zonotope.

### Input

- `M` -- square matrix with `size(M) == dim(P)`
- `P` -- sparse polynomial zonotope

### Output

The sparse polynomial zonotope resulting from applying the linear map.
"""
function linear_map(M::AbstractMatrix, P::SPZ)
    return SparsePolynomialZonotope(M * center(P),
                                    M * genmat_dep(P),
                                    M * genmat_indep(P),
                                    expmat(P),
                                    indexvector(P))
end
