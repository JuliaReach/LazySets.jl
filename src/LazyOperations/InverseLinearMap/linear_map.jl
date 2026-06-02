"""
    linear_map(M::AbstractMatrix, ilm::InverseLinearMap)

Return the linear map of a lazy inverse linear map.

### Input

- `M`   -- matrix
- `ilm` -- inverse linear map

### Output

The set representing the linear map of the lazy inverse linear map of a set.

### Notes

This implementation is inefficient because it computes the concrete inverse of
``M``, which is what `InverseLinearMap` is supposed to avoid.
"""
@validate function linear_map(M::AbstractMatrix, ilm::InverseLinearMap)
    return linear_map(M * inv(ilm.M), ilm.X)
end
