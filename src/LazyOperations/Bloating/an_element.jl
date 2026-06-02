"""
    an_element(B::Bloating)

Return some element of a bloated set.

### Input

- `B` -- bloated set

### Output

An element in the bloated set.

### Algorithm

This implementation disregards negative bloating and returns the result of
`an_element` for the wrapped set.
"""
function an_element(B::Bloating)
    if B.ε < 0
        throw(ArgumentError("negative bloating is not supported"))
    end
    return an_element(B.X)
end
