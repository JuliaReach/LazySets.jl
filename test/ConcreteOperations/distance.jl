for N in @tN([Float64, Float32, Rational{Int}])
    H1 = BallInf(N[0, 0], N(1))
    H2 = BallInf(N[-1 // 10, 0], N(1))
    H3 = Hyperrectangle(N[-16 // 10, -2], N[1 // 2, 1 // 2])
    S1 = Singleton(N[1, 1])
    S2 = Singleton(N[2, 3])

    for p in N[1, 2, Inf]
        @test distance(S1, S1; p=p) == N(0)
        @test distance(S1, S2; p=p) == distance(S2, S1; p=p) ≈ norm(N[1, 2], p)
        @test distance(H1, S1; p=p) == distance(S1, H1; p=p) == N(0)
        @test distance(H2, S1; p=p) == distance(S1, H2; p=p) ≈ norm(N[1 // 10, 0], p)
        @test distance(H3, S1; p=p) == distance(S1, H3; p=p) ≈ norm(N[21 // 10, 5 // 2], p)
        @test distance(H1, H2; p=p) == distance(H2, H1; p=p) == N(0)
        @test distance(H1, H3; p=p) == distance(H3, H1; p=p) ≈ norm(N[1 // 10, 1 // 2], p)
        @test distance(H2, H3; p=p) == distance(H3, H2; p=p) ≈ norm(N[0, 1 // 2], p)
    end
end
