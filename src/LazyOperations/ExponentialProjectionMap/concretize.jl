function concretize(epm::ExponentialProjectionMap)
    p = epm.projspmexp
    M = p.L * exp(Matrix(p.spmexp.M)) * p.R
    return linear_map(M, concretize(epm.X))
end
