"""
    ρ(d::AbstractVector, E::Ellipsoid)

Return the support function of an ellipsoid in a given direction.

### Input

- `d` -- direction
- `E` -- ellipsoid

### Output

The support function of the ellipsoid in the given direction.

### Algorithm

The support value is ``cᵀ d + ‖Bᵀ d‖₂``, where ``c`` is the center and
``Q = B Bᵀ`` is the shape matrix of `E`.
"""
function ρ(d::AbstractVector, E::Ellipsoid)
    return dot(E.center, d) + sqrt(inner(d, E.shape_matrix, d))
end
