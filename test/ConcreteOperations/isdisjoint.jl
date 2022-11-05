for N in [Float64, Float32, Rational{Int}]
    # using IA types
    X = IA.interval(N(0), N(1)) # IA
    Y = Interval(N(-1), N(2))
    Z = Interval(N(2), N(3))
    res, w = isdisjoint(X, Y, true)
    @test !isdisjoint(X, Y) && !isdisjoint(Y, X) && !res && w ∈ Interval(X) && w ∈ Y
    res, w = isdisjoint(X, Z, true)
    @test isdisjoint(X, Z) && isdisjoint(Z, X) && res && w == N[]

    X = IntervalBox(IA.interval(N(0), N(1)), IA.interval(N(0), N(1)))
    Y = Hyperrectangle(low=[N(-1), N(-1)], high=[N(2), N(2)])
    @test !isdisjoint(X, Y)
    @test !isdisjoint(Y, X)
end
