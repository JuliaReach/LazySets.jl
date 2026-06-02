"""
    an_element(rm::ResetMap)

Return some element of a reset map.

### Input

- `rm` -- reset map

### Output

An element in the reset map.

### Algorithm

This method relies on the `an_element` implementation for the wrapped set.
"""
function an_element(rm::ResetMap)
    return substitute(rm.resets, an_element(rm.X))
end
