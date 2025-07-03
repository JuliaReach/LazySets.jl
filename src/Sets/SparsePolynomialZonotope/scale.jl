function scale(α::Real, P::SparsePolynomialZonotope{N}) where {N<:Real}
    return SparsePolynomialZonotope(α * center(P),
                                    α * genmat_dep(P),
                                    α * genmat_indep(P),
                                    expmat(P),
                                    indexvector(P))
end