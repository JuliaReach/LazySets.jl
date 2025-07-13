@validate function affine_map(M::AbstractMatrix, E::Ellipsoid, v::AbstractVector)
    c = M * center(E)
    Q = _linear_map_shape_matrix(M, E)
    return Ellipsoid(c + v, Q)
end
