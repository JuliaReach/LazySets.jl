for N in @tN([Float64, Float32, Rational{Int}])
    # interval \ half-space
    x = Interval(N(1), N(3))
    U = HalfSpace(N[1], N(2))
    @test difference(x, U) == Interval(N(2), N(3))
    U = HalfSpace(N[-1], N(-2))
    @test difference(x, U) == Interval(N(1), N(2))
    for U in (HalfSpace(N[1], N(3)), HalfSpace(N[-1], N(-1)))
        @test difference(x, U) == EmptySet{N}(1)
    end
    for U in (HalfSpace(N[1], N(1)), HalfSpace(N[-1], N(-3)))
        @test difference(x, U) == x
    end
end
