function remove_small_generators(A::MatrixZonotope, tol::Real=1e-4; p::Real=Inf)
    Z = convert(Zonotope, A)
    Zred = remove_small_generators(Z, tol; p=p)

    # reshape to matrix zonotope
    gens = [Matrix(reshape(col, size(A))) for col in eachcol(genmat(Zred))]
    return MatrixZonotope(center(A), gens)
end
