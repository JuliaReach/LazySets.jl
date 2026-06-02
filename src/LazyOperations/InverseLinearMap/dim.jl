"""
    dim(ilm::InverseLinearMap)

Return the dimension of an inverse linear map.

### Input

- `ilm` -- inverse linear map

### Output

The ambient dimension of the inverse linear map.
"""
function dim(ilm::InverseLinearMap)
    return size(ilm.M, 1)
end
