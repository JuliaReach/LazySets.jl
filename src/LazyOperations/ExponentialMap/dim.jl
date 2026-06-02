"""
    dim(em::ExponentialMap)

Return the dimension of an exponential map.

### Input

- `em` -- exponential map

### Output

The ambient dimension of the exponential map.
"""
function dim(em::ExponentialMap)
    return size(em.expmat, 1)
end
