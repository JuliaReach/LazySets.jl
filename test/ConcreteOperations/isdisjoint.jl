for N in [Float64, Float32, Rational{Int}]

    # using IA types
    X = interval(N(0), N(1)) # IA
    Y = Interval(N(-1), N(2))
    @test !isdisjoint(X, Y)
    @test !isdisjoint(Y, X)

    X = IntevalBox(N(0) .. N(1), N(0) .. N(1))
    Y = Hyperrectangle(low=[N(-1), N(-1)], high=[N(2), N(2)])
    @test !isdisjoint(X, Y)
    @test !isdisjoint(Y, X)
end
