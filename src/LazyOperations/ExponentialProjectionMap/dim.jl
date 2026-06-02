"""
    dim(eprojmap::ExponentialProjectionMap)

Return the dimension of a projection of an exponential map.

### Input

- `eprojmap` -- projection of an exponential map

### Output

The ambient dimension of the projection of an exponential map.
"""
function dim(eprojmap::ExponentialProjectionMap)
    return size(eprojmap.projspmexp.L, 1)
end
