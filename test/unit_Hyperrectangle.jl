for N in [Float64, Rational{Int}, Float32]
    # 1D Hyperrectangle
    h = Hyperrectangle(N[0.], N[1.])
    # Test Dimension
    @test dim(h) == 1
    # Test Support Vector
    d = N[1.]
    @test σ(d, h) == N[1.]
    d = N[-1.]
    @test σ(d, h) == N[-1.]

    # 2D Hyperrectangle
    h = Hyperrectangle(N[0., 0.], N[1., 1.])
    # Test Dimension
    @test dim(h) == 2
    # Test Support Vector
    d = N[1., 1.]
    @test σ(d, h) == N[1., 1.]
    d = N[-1., 1.]
    @test σ(d, h) == N[-1., 1.]
    d = N[-1., -1.]
    @test σ(d, h) == N[-1., -1.]
    d = N[1., -1.]
    @test σ(d, h) == N[1., -1.]

    # 2D Hyperrectangle not 0-centered
    h = Hyperrectangle(N[1., 2.], N[1., 1.])
    # Test Dimension
    @test dim(h) == 2
    # Test Support Vector
    d = N[1., 1.]
    @test σ(d, h) == N[2., 3.]
    d = N[-1., 1.]
    @test σ(d, h) == N[0., 3.]
    d = N[-1., -1.]
    @test σ(d, h) == N[0., 1.]
    d = N[0., -1.]
    @test σ(d, h) == N[2., 1.]

    # 2D Hyperrectangle not same radius in each direction
    h = Hyperrectangle(N[0., 0.], N[1., 2.])
    # Test Dimension
    @test dim(h) == 2
    # Test Support Vector
    d = N[1., 1.]
    @test σ(d, h) == N[1., 2.]
    d = N[-1., 1.]
    @test σ(d, h) == N[-1., 2.]
    d = N[-1., -1.]
    @test σ(d, h) == N[-1., -2.]
    d = N[1., -1.]
    @test σ(d, h) == N[1., -2.]

    function isin(e, list)
    for x in list
        if x == e
        return true
        end
    end
    return false
    end

    # 2D Hyperrectangle not centered, not same radius, for vertex representation,
    # radius, and diameter
    h = Hyperrectangle(N[3., 2.], N[2., 1.])
    vl = vertices_list(h)
    # Test Vertices
    @test length(vl) == 4
    @test isin(N[1., 1.], vl) == true
    @test isin(N[1., 3.], vl) == true
    @test isin(N[5., 1.], vl) == true
    @test isin(N[5., 3.], vl) == true
    # norm
    @test norm(h) == norm(N[5., 3.], Inf)
    # radius
    @test radius(h) == norm(N[2., 1.], Inf)
    # diameter
    @test diameter(h) == norm(N[5., 3.] - N[1., 1.], Inf)

    # alternative constructors
    c = ones(N, 2)
    r = to_N(N, [0.1, 0.2])
    l = to_N(N, [0.9, 0.8])
    h = to_N(N, [1.1, 1.2])
    H1 = Hyperrectangle(c, r)
    H2 = Hyperrectangle(center=c, radius=r)
    H3 = Hyperrectangle(low=l, high=h)
    @test H1.center == H2.center
    @test H2.center ≈ H3.center
    @test H1.radius == H2.radius
    @test H2.radius ≈ H3.radius

    # Test low and high methods for a hyperrectangle
    H = Hyperrectangle(center=to_N(N, [-2.1, 5.6, 0.9]), radius=fill(to_N(N, 0.5), 3))
    @test low(H) ≈ to_N(N, [-2.6, 5.1, 0.4])
    @test high(H) ≈ to_N(N, [-1.6, 6.1, 1.4])

    # membership
    H = Hyperrectangle(N[1.0, 1.0], N[2.0, 3.0])
    @test !∈(N[-1.1, 4.1], H)
    @test ∈(N[-1.0, 4.0], H)
    @test !∈(N[-1.1, 4.1], H, N(0.09))
    @test ∈(N[-1.1, 4.1], H, N(0.11))

    # an_element function
    H = Hyperrectangle(N[1.0, 2.0], N[3.0, 4.0])
    @test an_element(H) ∈ H

    # subset
    H1 = Hyperrectangle(N[1.5, 1.5], N[0.5, 0.5])
    H2 = Hyperrectangle(N[2.0, 2.5], N[0.5, 0.5])
    H3 = Hyperrectangle(N[2.0, 2.0], N[2.0, 3.0])
    B1 = BallInf(N[2.0, 2.5], N(0.5))
    B2 = BallInf(N[2.0, 2.0], N(1.0))
    @test !⊆(H1, H2) && ⊆(H1, H3) && ⊆(H2, H3)
    subset, point = ⊆(H1, H2, true)
    @test !subset && point ∈ H1 && !(point ∈ H2)
    subset, point = ⊆(H1, H3, true)
    @test subset
    @test ⊆(H2, B1) && ⊆(B1, H2)
    @test ⊆(B1, B2) && !⊆(B2, B1)

    # intersection emptiness
    H1 = Hyperrectangle(N[1.0, 1.0], N[2.0, 2.0])
    H2 = Hyperrectangle(N[3.0, 3.0], N[2.0, 2.0])
    B1 = BallInf(N[2.0, 4.0], N(0.5))
    intersection_empty, point = is_intersection_empty(H1, H2, true)
    @test !is_intersection_empty(H1, H2) &&
    !intersection_empty && point ∈ H1 && point ∈ H2
    @test is_intersection_empty(H1, B1) && is_intersection_empty(H1, B1, true)[1]
end
