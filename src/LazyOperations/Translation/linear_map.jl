"""
    linear_map(M::AbstractMatrix, tr::Translation)

Concrete linear map of a translation.

### Input

- `M`  -- matrix
- `tr` -- translation of a set

### Output

A concrete set corresponding to the linear map.
The type of the result depends on the type of the set wrapped by `tr`.

### Algorithm

We compute `affine_map(M, tr.X, M * tr.v)`.
"""
@validate function linear_map(M::AbstractMatrix, tr::Translation)
    return affine_map(M, tr.X, M * tr.v)
end
