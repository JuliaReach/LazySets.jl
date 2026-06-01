"""
    center(cp::CartesianProduct)

Return the center of a Cartesian product of centrally-symmetric sets.

### Input

- `cp` -- Cartesian product of centrally-symmetric sets

### Output

The center of the Cartesian product.
"""
function center(cp::CartesianProduct)
    return vcat(center(cp.X), center(cp.Y))
end
