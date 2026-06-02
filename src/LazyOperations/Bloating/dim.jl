"""
    dim(B::Bloating)

Return the dimension of a bloated set.

### Input

- `B` -- bloated set

### Output

The ambient dimension of the bloated set.
"""
function dim(B::Bloating)
    return dim(B.X)
end
