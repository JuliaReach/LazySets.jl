function convert(::Type{DensePolynomialZonotope}, Z::AbstractZonotope{N}) where {N}
    c = center(Z)
    E = Matrix{N}[]
    F = Matrix{N}[]
    G = genmat(Z)
    return DensePolynomialZonotope(c, E, F, G)
end
