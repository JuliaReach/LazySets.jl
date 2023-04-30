for N in [Float64, Float32, Rational{Int}]
    # intervals, nonempty difference
    X = Interval(N(1), N(3))
    Y = Interval(N(2), N(3))
    @test minkowski_difference(X, Y) == Interval(N(-1), N(0))

    # intervals, empty difference
    X = Interval(N(1), N(3))
    Y = Interval(N(2), N(5))
    D = minkowski_difference(X, Y)
    @test D isa EmptySet{N} && dim(D) == 1

    # hyperrectangles, nonempty difference
    H1 = Hyperrectangle(N[2, 3], N[1, 3])
    H2 = Hyperrectangle(N[5/2, -1], N[1/2, 2])
    @test minkowski_difference(H1, H2) == Hyperrectangle(N[-1/2, 4], N[1/2, 1])

    # hyperrectangles, empty difference
    H1 = Hyperrectangle(N[2, 3], N[1, 3])
    H2 = Hyperrectangle(N[5, -1], N[2, 2])
    D = minkowski_difference(H1, H2)
    @test D isa EmptySet{N} && dim(D) == 2
end
