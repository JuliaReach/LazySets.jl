"""
# Extended help

    linear_map(M::AbstractMatrix, E::Ellipsoid)

### Algorithm

Given an ellipsoid ``⟨c, Q⟩`` and a matrix ``M``, the linear map yields the
ellipsoid ``⟨M c, M Q Mᵀ⟩``.
"""
function linear_map(M::AbstractMatrix, E::Ellipsoid)
    c = M * center(E)
    Q = _linear_map_shape_matrix(M, E)
    return Ellipsoid(c, Q)
end

function _linear_map_shape_matrix(M::AbstractMatrix, E::Ellipsoid)
    return M * shape_matrix(E) * M'
end
