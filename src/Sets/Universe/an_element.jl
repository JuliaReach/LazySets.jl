"""
    an_element(U::Universe{N}) where {N}

Return some element of a universe.

### Input

- `U` -- universe

### Output

The origin.
"""
function an_element(U::Universe{N}) where {N}
    return zeros(N, dim(U))
end
