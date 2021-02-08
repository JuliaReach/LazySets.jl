for N in [Float64, Float32]
    S = Star(N[3, 3], N[1 0; 0 1], BallInf(N[0, 0], N(1)))
    @test Star(N[3, 3], [N[1, 0], N[0, 1]], BallInf(N[0, 0], N(1))) == S
    @test center(S) == N[3, 3]
    @test basis(S) == N[1 0; 0 1]
    @test predicate(S) == BallInf(N[0, 0], N(1))
    @test dim(S) == 2
end
