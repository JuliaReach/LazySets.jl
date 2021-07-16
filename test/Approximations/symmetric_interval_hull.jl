for N in [Float64, Rational{Int}, Float32]
    S = Singleton(N[1, -2, 3, -4])
    H = Hyperrectangle(zeros(N, 4), N[1, 2, 3, 4])
    @test symmetric_interval_hull(S) == H

    M = Diagonal(ones(N, 4))
    @test symmetric_interval_hull(M * S) == H

    H2 = Hyperrectangle(zeros(N, 4), 2 * N[1, 2, 3, 4])
    @test symmetric_interval_hull(S + S) == H2
end
