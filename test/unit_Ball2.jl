for N in [Float64, Float32]
    # 1D Ball2
    b = Ball2(N[0.], N(1.))
    # Test Dimension
    @test dim(b) == 1
    # Test Support Vector
    d = N[1.]
    @test σ(d, b) == N[1.]
    d = N[-1.]
    @test σ(d, b) == N[-1.]

    # 2D Ball2
    b = Ball2(N[0., 0.], N(1.))
    # Test Dimension
    @test dim(b) == 2
    # Test Support Vector
    d = N[1., 0.]
    @test σ(d, b) == N[1., 0.]
    d = N[-1., 0.]
    @test σ(d, b) == N[-1., 0.]
    d = N[0., 1.]
    @test σ(d, b) == N[0., 1.]
    d = N[0., -1.]
    @test σ(d, b) == N[0., -1.]
    d = N[0., 0.]
    @test σ(d, b) == N[0., 0.]

    # 2D Ball2 not 0-centered
    b = Ball2(N[1., 2.], N(1.))
    # Test Dimension
    @test dim(b) == 2
    # Test Support Vector
    d = N[1., 0.]
    @test σ(d, b) == N[2., 2.]
    d = N[-1., 0.]
    @test σ(d, b) == N[0., 2.]
    d = N[0., 1.]
    @test σ(d, b) == N[1., 3.]
    d = N[0., -1.]
    @test σ(d, b) == N[1., 1.]

    # 2D Ball2 radius =/= 1
    b = Ball2(N[0., 0.], N(2.))
    # Test Dimension
    @test dim(b) == 2
    # Test Support Vector
    d = N[1., 0.]
    @test σ(d, b) == N[2., 0.]
    d = N[-1., 0.]
    @test σ(d, b) == N[-2., 0.]
    d = N[0., 1.]
    @test σ(d, b) == N[0., 2.]
    d = N[0., -1.]
    @test σ(d, b) == N[0., -2.]

    # an_element function
    b = Ball2(N[1., 2.], N(2.))
    @test an_element(b) ∈ b

    # membership
    x = N[0., 0.]
    b = Ball2(N[1., 1.], N(0.9))
    @test !in(x, b)
    @test !in(x, b, N(0.5))
    @test in(x, b, N(0.55))

    # subset
    b1 = Ball2(N[1., 2.], N(2.))
    b2 = Ball2(N[1., 2.], N(0.))
    b3 = Ball2(N[1.7, 2.7], N(1.))
    s = Singleton(N[1., 2.])
    subset, point = ⊆(b1, s, true)
    @test !⊆(b1, s) && !subset && point ∈ b1 && !(point ∈ s)
    @test ⊆(b2, s) && ⊆(b2, s, true)[1]
    @test !⊆(b1, BallInf(N[1., 2.], N(1.)))
    @test ⊆(b2, BallInf(N[1., 2.], N(2.)))
    subset, point = ⊆(b1, b3, true)
    @test !⊆(b1, b3) && !subset && point ∈ b1 && !(point ∈ b3)
    @test ⊆(b3, b1) && ⊆(b3, b1, true)[1]

    # intersection
    b1 = Ball2(N[0., 0.], N(2.))
    b2 = Ball2(N[2., 2.], N(2.))
    b3 = Ball2(N[4., 4.], N(2.))
    intersection_empty, point = is_intersection_empty(b1, b2, true)
    @test !is_intersection_empty(b1, b2) && !intersection_empty && point ∈ b1 && point ∈ b2
    @test is_intersection_empty(b1, b3) && is_intersection_empty(b1, b3, true)[1]
end
