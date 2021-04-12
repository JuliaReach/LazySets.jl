for N in [Float64, Float32, Rational{Int}]

    # constructor with basis matrix
    S = Star(N[3, 3], N[1 0; 0 1], BallInf(N[0, 0], N(1)))

    # constructor vector-of-vectors basis
    W = Star(N[3, 3], [N[1, 0], N[0, 1]], BallInf(N[0, 0], N(1)))

    @test S ≈ W

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
    @test σ(N[-1, -1], S) == N[2, 2]
    @test σ(N[1, -1], S) == N[4, 2]
    @test σ(N[-1, 1], S) == N[2, 4]

    # conversion from a polyhedral set
    P = HPolygon([HalfSpace(N[1, 1], N(-0.5)),
                  HalfSpace(N[-1.5, 0.5], N(4)),
                  HalfSpace(N[-0.5, -0.05], N(1.5)),
                  HalfSpace(N[1.5, -0.25], N(-1.5))])
    S = convert(STAR, P)
    @test basis(S) == N[1 0; 0 1]
    @test center(S) == N[0, 0]
    @test predicate(S) == P

    # intersection with HalfSpace
    H = HalfSpace(N[0, 1], N(0))
    I = intersection(S, H)
    addconstraint!(P, H)
    @test isequivalent(I, P)
end
