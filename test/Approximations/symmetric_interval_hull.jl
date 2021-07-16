for N in [Float64, Rational{Int}, Float32]
    # singleton
    S = Singleton(N[1, -2, 3, -4])
    H = Hyperrectangle(zeros(N, 4), N[1, 2, 3, 4])
    @test symmetric_interval_hull(S) == H

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
end
