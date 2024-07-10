"""
    dim(U::Universe)

Return the dimension of a universe.

### Input

- `U` -- universe

### Output

The ambient dimension of a universe.
"""
function dim(U::Universe)
    return U.dim
end
