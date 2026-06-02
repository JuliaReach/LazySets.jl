function concretize(am::AffineMap)
    return affine_map(am.M, concretize(am.X), am.v)
end
