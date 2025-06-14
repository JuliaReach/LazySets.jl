function linear_map(M::AbstractMatrix, P::SimpleSparsePolynomialZonotope)
    return SimpleSparsePolynomialZonotope(M * center(P), M * genmat(P), expmat(P))
end
