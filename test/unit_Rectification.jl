for N in [Float64, Rational{Int}, Float32]
    I1 = Interval(N(-1), N(1))
    I2 = Interval(N(2), N(3))
    I3 = Interval(N(-2), N(-1))
    B1 = BallInf(N[1, 1], N(2))
    B2 = BallInf(N[3, 1], N(2))
    B3 = BallInf(N[3, 3], N(2))

    # constructor
    RI1 = Rectification(I1)
    RI2 = Rectification(I2)
    RI3 = Rectification(I3)
    RB1 = Rectification(B1)
    RB2 = Rectification(B2)
    RB3 = Rectification(B3)

    # dimension
    @test dim(RI1) == dim(RI2) == dim(RI3) == 1
    @test dim(RB1) == dim(RB2) == dim(RB3) == 2

    # support vector
    # hyperrectangles
    @test σ(N[1], RI1) == N[1] && σ(N[-1], RI1) == N[0]
    @test σ(N[1], RI2) == N[3] && σ(N[-1], RI2) == N[2]
    @test σ(N[1], RI3) == N[0] && σ(N[-1], RI3) == N[0]
    @test σ(N[1, 1], RB1) == N[3, 3] && σ(N[-1, 1], RB1) == N[0, 3]
    @test σ(N[1, 1], RB2) == N[5, 3] && σ(N[-1, 1], RB2) == N[1, 3]
    @test σ(N[1, 1], RB3) == N[5, 5] && σ(N[-1, 1], RB3) == N[1, 5]
    # other sets in 1D fall back to interval conversion
    @test σ(N[1], Rectification(Ball1(N[0], N(1)))) == N[1]
    # other sets in higher dimensions throw an error
    @test_throws ErrorException σ(N[1, 1], Rectification(Ball1(N[0, 0], N(1))))

    # an_element
    x = an_element(RI1)
    @test x isa Vector{N} && length(x) == 1 && N(0) <= x[1] <= N(1)
    x = an_element(RB2)
    @test x isa Vector{N} && length(x) == 2 && N(1) <= x[1] <= N(3) &&
          N(0) <= x[2] <= N(3)

    # membership
    @test N[-1, 1] ∉ RB1
    @test N[1, 1] ∈ RB1
    @test N[0, 1] ∈ RB1
    @test_throws ErrorException N[0, 4] ∈ RB1

    # emptiness
    @test !isempty(RI1)
    @test isempty(Rectification(EmptySet{N}()))
end
