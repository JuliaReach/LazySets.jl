"""
    volume(cp::CartesianProduct)

Compute the volume of a Cartesian product.

### Input

- `cp` -- Cartesian product

### Output

The volume.
"""
function volume(cp::CartesianProduct)
    return volume(cp.X) * volume(cp.Y)
end
