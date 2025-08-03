"""
    expmat(P::SimpleSparsePolynomialZonotope)

Return the matrix of exponents of a simple sparse polynomial zonotope.

### Input

- `P` -- simple sparse polynomial zonotope

### Output

The matrix of exponents, where each column is a multidegree.

### Notes

In the exponent matrix, each row corresponds to a parameter (``\alpha_k`` in the
mathematical set definition) and each column corresponds to a monomial.

### Examples

```jldoctest
julia> S = SimpleSparsePolynomialZonotope([2.0, 0], [1 2;2 2.], [1 4;1 2])
SimpleSparsePolynomialZonotope{Float64, Vector{Float64}, Matrix{Float64}, Matrix{Int64}}([2.0, 0.0], [1.0 2.0; 2.0 2.0], [1 4; 1 2])

julia> expmat(S)
2×2 Matrix{Int64}:
 1  4
 1  2
```
"""
function expmat(P::SimpleSparsePolynomialZonotope)
    return P.E
end
