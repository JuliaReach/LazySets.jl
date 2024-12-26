function linear_map(M::AbstractMatrix, P::SPZ)
    return SparsePolynomialZonotope(M * center(P),
                                    M * genmat_dep(P),
                                    M * genmat_indep(P),
                                    expmat(P),
                                    indexvector(P))
end
