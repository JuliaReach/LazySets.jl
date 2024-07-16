"""
    affine_map(M::AbstractMatrix, X::Star, v::AbstractVector)

Return the affine map of a star.

### Input

- `M` -- matrix
- `X` -- star
- `v` -- vector

### Output

The star obtained by applying the affine map with matrix `M` and displacement
`v` to `X`.
"""
function affine_map(M::AbstractMatrix, X::Star, v::AbstractVector)
    c′ = M * X.c + v
    V′ = M * X.V
    P′ = X.P
    return Star(c′, V′, P′)
end
