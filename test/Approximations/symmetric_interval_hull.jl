for N in [Float64, Rational{Int}, Float32]
    # singleton
    S = Singleton(N[1, -2, 3, -4])
    H = Hyperrectangle(zeros(N, 4), N[1, 2, 3, 4])
    @test symmetric_interval_hull(S) == H
    @test box_approximation_symmetric(S) == H

    # linear map of singleton
    M = Diagonal(ones(N, 4))
    @test symmetric_interval_hull(M * S) == H

    # Minkowski sum of singletons
    H2 = Hyperrectangle(zeros(N, 4), 2 * N[1, 2, 3, 4])
    @test symmetric_interval_hull(S + S) == H2

    # linear map of hyperrectangle
    M = Diagonal(N[1, -1, -2, 2])
    H1 = Hyperrectangle(5 * ones(N, 4), N[1, 2, 3, 4])
    H2 = Hyperrectangle(zeros(N, 4), N[6, 7, 16, 18])
    @test symmetric_interval_hull(M * H1) == H2

    # LineSegment & V-rep
    v = [N[1, 0], [0, 2]]
    for V in (LineSegment(v[1], v[2]), VPolygon(v), VPolytope(v))
        @test symmetric_interval_hull(V) == Hyperrectangle(zeros(N, 2), N[1, 2])
    end
end

# tests that only work with Float64 and Float32
for N in [Float64, Float32]
    # exponential map of singleton
    M = sparse(diagm(N[1, 1, 1]))
    E = SparseMatrixExp(M) * Singleton(N[1, 2, 3])
    H = Hyperrectangle(zeros(N, 3), N[ℯ, 2ℯ, 3ℯ])
    @test symmetric_interval_hull(E) ≈ H

    # exponential map of hyperrectangle
    M = sparse(diagm(N[1, 1, 1]))
    E = SparseMatrixExp(M) * Hyperrectangle(N[1, 2, 3], N[1, 2, 3])
    H = symmetric_interval_hull(E)
    @test H ≈ Hyperrectangle(zeros(N, 3), N[2ℯ, 4ℯ, 6ℯ])
    # matrix with negative entries
    M = sparse(diagm(N[1, -1]))
    E = SparseMatrixExp(M) * Hyperrectangle(N[1, 1], N[1, 1])
    H = symmetric_interval_hull(E)
    @test H ≈ Hyperrectangle(zeros(N, 2), N[2ℯ, 2ℯ^(-1)])
end
