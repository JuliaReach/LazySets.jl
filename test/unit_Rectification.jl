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

    # set
    @test set(RI1) == I1
    @test set(RI2) == I2
    @test set(RI3) == I3
    @test set(RB1) == B1
    @test set(RB2) == B2
    @test set(RB3) == B3

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
    # Cartesian products
    RC1 = Rectification(I1 × I2)
    @test σ(N[1, 1], RC1) == N[1, 3] && σ(N[-1, -1], RC1) == N[0, 2]
    RC2 = Rectification(CartesianProductArray([I1, I2, I3]))
    @test σ(N[1, 1, 1], RC2) == N[1, 3, 0] &&
          σ(N[-1, -1, -1], RC2) == N[0, 2, 0]
    # other sets in 1D fall back to interval conversion
    @test σ(N[1], Rectification(Ball1(N[0], N(1)))) == N[1]
    # other sets in higher dimensions compute lazy linear maps and intersections
    @test σ(N[-1, -1], Rectification(Ball1(N[-2, -2], N(1)))) == N[0, 0]
    # this may throw an error due to a lazy intersection
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
    @test isempty(Rectification(EmptySet{N}(1)))

    # boundedness
    @test isbounded(RI1)
    @test isbounded(Rectification(HalfSpace(N[1], N(0))))
    @test !isbounded(Rectification(Universe{N}(2)))
    @test !isbounded(Rectification(HalfSpace(N[-1], N(0))))

    # conversion
    @test convert(Interval, RI1) == Interval(N(0), N(1))
    @test convert(Hyperrectangle, RB1) ==
                  Hyperrectangle(low=N[0, 0], high=N[3, 3])

    # support function
    # upper right quadrant
    B = Ball1(N[2, 2], N(1))
    RB = Rectification(B)
    for d in [N[1, 0], N[-1, 0], N[0, 1], N[0, -1]]
        @test ρ(d, RB) == ρ(d, B)
    end
    # lower right quadrant
    B = Ball1(N[2, -2], N(1))
    RB = Rectification(B)
    for d in [N[1, 0], N[-1, 0]]
        @test ρ(d, RB) == ρ(d, B)
    end
    @test ρ(N[0, 1], RB) == ρ(N[0, -1], RB) == N(0)
    # upper left quadrant
    B = Ball1(N[-2, 2], N(1))
    RB = Rectification(B)
    for d in [N[0, 1], N[0, -1]]
        @test ρ(d, RB) == ρ(d, B)
    end
    @test ρ(N[1, 0], RB) == ρ(N[-1, 0], RB) == N(0)
    # lower left quadrant
    B = Ball1(N[-2, -2], N(1))
    RB = Rectification(B)
    @test ρ(N[1, 0], RB) == ρ(N[-1, 0], RB) == ρ(N[0, 1], RB) ==
          ρ(N[0, -1], RB) == N(0)

    # concretize
    @test concretize(RI1) == rectify(I1)
    @test concretize(RI2) == rectify(I2)
    @test concretize(RI3) == rectify(I3)
end

# tests that only work with Float64
for N in [Float64]
    # support function
    # some of the tests do not work because of insufficient precision in the
    # intersection; if the precision changes, these tests can be replaced by
    # their true (commented out) expected results

    # two right quadrants
    B = Ball1(N[2, 0], N(1))
    RB = Rectification(B)
    for d in [N[1, 0], N[-1, 0], N[0, 1]]
        @test ρ(d, RB) ≈ ρ(d, B)
    end
#     @test ρ(N[0, -1], RB) ≈ N(0)
    @test N(0) ≤ ρ(N[0, -1], RB) ≤ N(1e-9)
    # all four quadrants
    P = VPolygon([N[-1, 1], N[-1.5, 0.5], N[1.5, 0.5], N[1, -0.5]])
    RP = Rectification(P)
    @test ρ(N[1, 0], RP) ≈ N(1.5)
    @test ρ(N[0, 1], RP) ≈ N(1)
    @test ρ(N[1, 1], RP) ≈ ρ(N[1, 1], P) == N(2)
#     @test ρ(N[-1, 0], RP) ≈ N(0)
    @test N(0) ≤ ρ(N[-1, 0], RP) ≤ N(1e-8)
#     @test ρ(N[0, -1], RP) ≈ N(0)
    @test N(0) ≤ ρ(N[0, -1], RP) ≤ N(1e-8)
#     @test ρ(N[-1, -1], RP) ≈ N(0)
    @test N(0) ≤ ρ(N[-1, -1], RP) ≤ N(1.1e-1)
end
