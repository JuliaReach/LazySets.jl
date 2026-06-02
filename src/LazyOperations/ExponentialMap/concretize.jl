function concretize(em::ExponentialMap)
    return exponential_map(Matrix(em.expmat.M), concretize(em.X))
end
