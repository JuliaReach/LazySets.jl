function minkowski_sum(A::MatrixZonotope, B::MatrixZonotope) 
    c = center(A) + center(B)
    gens = vcat(generators(A), generators(B))
    return MatrixZonotope(c, gens)
end
