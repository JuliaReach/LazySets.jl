"""
    dim(lm::LinearMap)

Return the dimension of a linear map.

### Input

- `lm` -- linear map

### Output

The ambient dimension of the linear map.
"""
function dim(lm::LinearMap)
    return size(lm.M, 1)
end
