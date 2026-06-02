"""
    linear_map(M::AbstractMatrix, lm::LinearMap)

Return the linear map of a lazy linear map.

### Input

- `M`  -- matrix
- `lm` -- linear map

### Output

A set representing the linear map.
"""
@validate function linear_map(M::AbstractMatrix, lm::LinearMap)
    return linear_map(M * lm.M, lm.X)
end
