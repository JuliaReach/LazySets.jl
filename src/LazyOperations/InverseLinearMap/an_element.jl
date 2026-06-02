"""
    an_element(ilm::InverseLinearMap)

Return some element of an inverse linear map.

### Input

- `ilm` -- inverse linear map

### Output

An element in the inverse linear map.
It relies on the `an_element` function of the wrapped set.
"""
function an_element(lm::InverseLinearMap)
    return lm.M \ an_element(lm.X)
end
