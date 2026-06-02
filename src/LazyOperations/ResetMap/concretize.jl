function concretize(rm::ResetMap)
    return affine_map(matrix(rm), concretize(set(rm)), vector(rm))
end
