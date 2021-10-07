for N in [Float64, Float32, Rational{Int}]

    # using IA types
    X = IA.interval(N(0), N(1)) # IA
    Y = IA.Interval(N(-1), N(2))
    @test !isdisjoint(X, Y)
    @test !isdisjoint(Y, X)

    X = IntervalBox(IA.interval(N(0), N(1)), IA.interval(N(0), N(1)))
    Y = Hyperrectangle(low=[N(-1), N(-1)], high=[N(2), N(2)])
    @test !isdisjoint(X, Y)
    @test !isdisjoint(Y, X)
end
