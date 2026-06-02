"""
    an_element(lm::LinearMap)

Return some element of a linear map.

### Input

- `lm` -- linear map

### Output

An element in the linear map.
It relies on the `an_element` function of the wrapped set.
"""
function an_element(lm::LinearMap)
    return lm.M * an_element(lm.X)
end
