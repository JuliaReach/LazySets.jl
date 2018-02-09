for N in [Float64, Rational{Int}, Float32]
    # Sum of 2D centered balls in norm 1 and infinity
    b1 = BallInf(N[0., 0.], N(2.))
    b2 = Ball1(N[0., 0.], N(1.))
    # Test Construction
    X = MinkowskiSum(b1, b2)
    @test X.X == b1
    @test X.Y == b2
    # Test Dimension
    @test dim(X) == 2
    # Test Support Vector
    d = N[1., 0.]
    v = σ(d, X)
    @test v[1] == N(3.)
    d = N[-1., 0.]
    v = σ(d, X)
    @test v[1] == N(-3.)
    d = N[0., 1.]
    v = σ(d, X)
    @test v[2] == N(3.)
    d = N[0., -1.]
    v = σ(d, X)
    @test v[2] == N(-3.)

    # Sum of not-centered 2D balls in norm 1 and infinity
    b1 = BallInf(N[-1., 3.], N(2.))
    b2 = Ball1(N[1., 2.], N(1.))
    s = b1 + b2
    # Test Support Vector
    d = N[1., 0.]
    v = σ(d, s)
    @test v[1] == N(3.)
    d = N[-1., 0.]
    v = σ(d, s)
    @test v[1] == N(-3.)
    d = N[0., 1.]
    v = σ(d, s)
    @test v[2] == N(8.)
    d = N[0., -1.]
    v = σ(d, s)
    @test v[2] == N(2.)

    # Sum of array of LazySet
    # 2-elements
    ms = MinkowskiSum(Singleton(N[1.]), Singleton(N[2.]))
    @test ρ(N[1.], ms) == 3.
    @test ρ(N[-1.], ms) == -3.
    # 3-elements
    ms = MinkowskiSum(Singleton(N[1.]), MinkowskiSum(Singleton(N[2.]), Singleton(N[3.])))
    @test ρ(N[1.], ms) == N(6.)
    @test ρ(N[-1.], ms) == N(-6.)

    # array getter
    v = Vector{N}(0)
    @test array(MinkowskiSumArray()) == v

    # neutral and absorbing element
    z = ZeroSet{N}(2)
    e = EmptySet{N}()
    b = BallInf(N[0., 0.], N(2.))
    @test b + z == z + b == b
    @test b + e == e + b == e + e == e
    @test z + e == e + z == e
end
