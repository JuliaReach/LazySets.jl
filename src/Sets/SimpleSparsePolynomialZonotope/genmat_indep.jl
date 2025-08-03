function genmat_indep(P::SimpleSparsePolynomialZonotope{N}) where {N}
    return Matrix{N}(undef, dim(P), 0)
end
