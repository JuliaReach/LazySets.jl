"""
    dim(Z::ZeroSet)

Return the ambient dimension of a zero set.

### Input

- `Z` -- zero set

### Output

The ambient dimension of the zero set.
"""
function dim(Z::ZeroSet)
    return Z.dim
end
