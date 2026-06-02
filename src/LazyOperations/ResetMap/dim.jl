"""
    dim(rm::ResetMap)

Return the dimension of a reset map.

### Input

- `rm` -- reset map

### Output

The ambient dimension of a reset map.
"""
function dim(rm::ResetMap)
    return dim(rm.X)
end
