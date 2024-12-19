"""
# Extended help

    ρ(d::AbstractVector, E::Ellipsoid)

### Algorithm

The support value is ``cᵀ d + ‖Bᵀ d‖₂``, where ``c`` is the center and
``Q = B Bᵀ`` is the shape matrix of `E`.
"""
function ρ(d::AbstractVector, E::Ellipsoid)
    return dot(E.center, d) + sqrt(inner(d, E.shape_matrix, d))
end
