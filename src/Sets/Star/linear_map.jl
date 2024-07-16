"""
    linear_map(M::AbstractMatrix, X::Star)

Return the linear map of a star.

### Input

- `M` -- matrix
- `X` -- star

### Output

The star obtained by applying `M` to `X`.
"""
function linear_map(M::AbstractMatrix, X::Star)
    c′ = M * X.c
    V′ = M * X.V
    P′ = X.P
    return Star(c′, V′, P′)
end
