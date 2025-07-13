@validate function linear_map(M::AbstractMatrix, P::DensePolynomialZonotope)
    c = M * P.c
    E = [M * Ei for Ei in P.E]
    F = [M * Fi for Fi in P.F]
    G = M * P.G
    return DensePolynomialZonotope(c, E, F, G)
end
