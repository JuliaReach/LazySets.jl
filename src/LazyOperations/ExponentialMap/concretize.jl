function concretize(em::ExponentialMap)
    @assert em isa (ExponentialMap{N,S,<:SparseMatrixExp} where {N,S}) "matrix sets are not supported"
    return exponential_map(Matrix(em.expmat.M), concretize(em.X))
end
