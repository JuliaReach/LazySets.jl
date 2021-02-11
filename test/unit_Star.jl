for N in [Float64, Float32, Rational{Int}]

    # constructor with basis matrix
    S = Star(N[3, 3], N[1 0; 0 1], BallInf(N[0, 0], N(1)))

    # constructor vector-of-vectors basis
    W = Star(N[3, 3], [N[1, 0], N[0, 1]], BallInf(N[0, 0], N(1)))

    # S == W FIXME requires Star <: LazySet
    @test S.c == W.c && S.V == W.V && S.P == W.P

    # getter functions
    @test center(S) == N[3, 3]
    @test basis(S) == N[1 0; 0 1]
    @test predicate(S) == BallInf(N[0, 0], N(1))
    @test dim(S) == 2

    # support function
    @test ρ(N[1, 0], S) == N(4)
    @test -ρ(N[-1, 0], S) == N(2)
    @test ρ(N[0, 1], S) == N(4)
    @test -ρ(N[0, -1], S) == N(2)

    # support vector
    @test σ(N[1, 1], S) == N[4, 4]
    @test -σ(N[-1, -1], S) == N[-2, 2]
    @test σ(N[1, -1], S) == N[4, 2]
    @test σ(N[-1, 1], S) == N[2, 4]
end
