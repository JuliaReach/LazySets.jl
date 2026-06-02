"""
    dim(cp::CartesianProduct)

Return the dimension of a Cartesian product.

### Input

- `cp` -- Cartesian product

### Output

The ambient dimension of the Cartesian product.
"""
function dim(cp::CartesianProduct)
    return dim(cp.X) + dim(cp.Y)
end
