"""
    dim(qm::QuadraticMap)

Return the dimension of a quadratic map.

### Input

- `qm` -- quadratic map

### Output

The ambient dimension of the quadratic map.
"""
function dim(qm::QuadraticMap)
    return dim(qm.X)
end
