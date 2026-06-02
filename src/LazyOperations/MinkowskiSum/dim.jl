"""
    dim(ms::MinkowskiSum)

Return the dimension of a Minkowski sum of two sets.

### Input

- `ms` -- Minkowski sum of two sets

### Output

The ambient dimension of the Minkowski sum of two sets.
"""
function dim(ms::MinkowskiSum)
    return dim(ms.X)
end
