function polynomial_order(P::SimpleSparsePolynomialZonotope)
    return maximum(sum, eachcol(expmat(P)))
end
