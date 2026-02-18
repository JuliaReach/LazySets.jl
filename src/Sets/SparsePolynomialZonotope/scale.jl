function scale!(α::Real, P::SparsePolynomialZonotope)
    P.c .*= α
    P.G .*= α
    P.GI .*= α
    return P
end

function scale(α::Real, P::SparsePolynomialZonotope)
    return scale!(α, copy(P))
end
