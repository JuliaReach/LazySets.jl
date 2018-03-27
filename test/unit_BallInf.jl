for N in [Float64, Rational{Int}, Float32]
    # 1D BallInf
    b = BallInf(N[0.], N(1.))
    # Test Dimension
    @test dim(b) == 1
    # Test Support Vector
    d = N[1.]
    @test σ(d, b) == N[1.]
    d = N[-1.]
    @test σ(d, b) == N[-1.]

    # 2D BallInf
    b = BallInf(N[0., 0.], N(1.))
    # Test Dimension
    @test dim(b) == 2
    # Test Support Vector
    d = N[1., 1.]
    @test σ(d, b) == N[1., 1.]
    d = N[-1., 1.]
    @test σ(d, b) == N[-1., 1.]
    d = N[-1., -1.]
    @test σ(d, b) == N[-1., -1.]
    d = N[1., -1.]
    @test σ(d, b) == N[1., -1.]

    # 2D BallInf not 0-centered
    b = BallInf(N[1., 2.], N(1.))
    # Test Dimension
    @test dim(b) == 2
    # Test Support Vector
    d = N[1., 1.]
    @test σ(d, b) == N[2., 3.]
    d = N[-1., 1.]
    @test σ(d, b) == N[0., 3.]
    d = N[-1., -1.]
    @test σ(d, b) == N[0., 1.]
    d = N[0., -1.]
    @test σ(d, b) == N[2., 1.]

    # 2D BallInf radius =/= 1
    b = BallInf(N[0., 0.], N(2.))
    # Test Dimension
    @test dim(b) == 2
    # Test Support Vector
    d = N[1., 1.]
    @test σ(d, b) == N[2., 2.]
    d = N[-1., 1.]
    @test σ(d, b) == N[-2., 2.]
    d = N[-1., -1.]
    @test σ(d, b) == N[-2., -2.]
    d = N[1., -1.]
    @test σ(d, b) == N[2., -2.]

    # membership
    b = BallInf(N[1., 1.], N(1.))
    @test !∈(N[.5, -.5], b)
    @test ∈(N[.5, 1.5], b)
    @test !∈(N[.5, -.5], b, N(0.4))
    @test ∈(N[.5, -.5], b, N(0.6))

    # an_element function
    b = BallInf(N[1., 2.], N(3.))
    @test an_element(b) ∈ b
end
