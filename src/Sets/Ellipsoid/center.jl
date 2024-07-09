"""
    center(E::Ellipsoid)

Return the center of the ellipsoid.

### Input

- `E` -- ellipsoid

### Output

The center of the ellipsoid.
"""
function center(E::Ellipsoid)
    return E.center
end
