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
end
