function affine_map(M::AbstractMatrix, E::Ellipsoid, v::AbstractVector)
    c = _linear_map_center(M, E)
    Q = _linear_map_shape_matrix(M, E)
    return Ellipsoid(c + v, Q)
end
