function remove_redundant_generators(A::MatrixZonotope)
    Z = convert(Zonotope, A)
    Zred = remove_redundant_generators(Z)

    # reshape to matrix zonotope
    gens = [Matrix(reshape(col, size(A))) for col in eachcol(genmat(Zred))]
    return MatrixZonotope(center(A), gens)
end
