function scale(α::Real, P::DensePolynomialZonotope)
    return DensePolynomialZonotope(α * P.c, α * P.E, α * P.F, α * P.G)
end

function scale!(α::Real, P::DensePolynomialZonotope)
    P.c .*= α
    P.E .*= α
    P.F .*= α
    P.G .*= α
    return P
end
