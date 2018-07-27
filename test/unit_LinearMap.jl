for N in [Float64, Rational{Int}, Float32]
    # π/2 trigonometric rotation
    b = BallInf(N[1., 2.], N(1.))
    M = N[0. -1. ; 1. 0.]
    # Test Construction
    lm1 = LinearMap(M, b)
    @test lm1.M == M
    @test lm1.X == b
    # Test Dimension
    @test dim(lm1) == 2
    # Test Support Vector
    d = N[1., 1.]
    @test σ(d, lm1) == N[-1., 2.]
    d = N[-1., 1.]
    @test σ(d, lm1) == N[-3., 2.]
    d = N[-1., -1.]
    @test σ(d, lm1) == N[-3., 0.]
    d = N[1., -1.]
    @test σ(d, lm1) == N[-1., 0.]

    # 2D -> 1D Projection
    b = BallInf(N[1., 2.], N(1.))
    M = N[1. 0.]
    lm = M*b
    # Test Dimension
    @test dim(lm) == 1
    # Test Support Vector
    d = N[1.]
    @test σ(d, lm) == N[2.]
    d = N[-1.]
    @test σ(d, lm) == N[0.]

    # scalar multiplication
    b = BallInf(N[0., 0.], N(1.))
    lm = N(2.) * b
    # repeated scalar multiplication
    lm2 = N(2.) * lm
    # Test Dimension
    @test dim(lm) == 2
    # Test Support Vector
    d = N[1., 1.]
    @test σ(d, lm) == N[2., 2.]
    d = N[-1., 1.]
    @test σ(d, lm) == N[-2., 2.]
    d = N[-1., -1.]
    @test σ(d, lm) == N[-2., -2.]
    d = N[1., -1.]
    @test σ(d, lm) == N[2., -2.]

    # Nested construction
    lm1_copy = LinearMap(eye(N, 2), lm1)
    @test lm1_copy.M == lm1.M
    @test lm1_copy.X == lm1.X

    # an_element function
    lm = N(2.) * BallInf(N[0., 0.], N(1.))
    an_element(lm)
    @test an_element(lm) ∈ lm

    # check linear map between vector and set
    X = Interval([0.9, 1.10445])
    a = [-1., 2.]
    @test a * X isa LinearMap
    X = BallInf(N[1.0], N(1.10445))
    a = N[-1, 2.]
    @test a * X isa LinearMap{N, N}

    # linear map with a ZeroSet
    X = N[0. -1. ; 1. 0.] * ZeroSet{N}(2)
    @test X isa ZeroSet{N} && dim(X) == 2
end
