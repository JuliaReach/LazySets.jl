"""
    center(ms::MinkowskiSum)

Return the center of a Minkowski sum of two centrally-symmetric sets.

### Input

- `ms` -- Minkowski sum of two centrally-symmetric sets

### Output

The center of the Minkowski sum.
"""
function center(ms::MinkowskiSum)
    return center(ms.X) + center(ms.Y)
end
