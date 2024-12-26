function translate(S::SparsePolynomialZonotope, v::AbstractVector)
    c = center(S) + v
    return SparsePolynomialZonotope(c, genmat_dep(S), genmat_indep(S), expmat(S))
end
