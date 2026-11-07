@doc """
    affine_map(M::AbstractMatrix, X::LazySet, v::AbstractVector)

Compute the affine map ``M · X + v``.

### Input

- `M` -- matrix
- `X` -- set
- `v` -- translation vector

### Output

A set representing the affine map ``M · X + v``.
"""
affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector)
function affine_map end
