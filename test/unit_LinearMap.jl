for N in [Float64, Rational{Int}, Float32]
    # π/2 trigonometric rotation
    b = BallInf(N[1., 2.], N(1.))
    M = N[0. -1. ; 1. 0.]
    # Test Construction
    lm1 = LinearMap(M, b)
    @test lm1.M == M
    @test lm1.sf == b
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
    @test lm1_copy.sf == lm1.sf
end
