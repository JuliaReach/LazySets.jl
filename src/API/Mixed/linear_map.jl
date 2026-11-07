@doc """
    linear_map(M::AbstractMatrix, X::LazySet)

Compute the linear map ``M · X``.

### Input

- `M` -- matrix
- `X` -- set

### Output

A set representing the linear map ``M · X``.
"""
linear_map(::AbstractMatrix, ::LazySet)
function linear_map end
