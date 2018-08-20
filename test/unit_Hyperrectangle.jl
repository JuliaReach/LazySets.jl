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

    # 2D Hyperrectangle not centered, not same radius, for vertex representation,
    # radius, and diameter
    h = Hyperrectangle(N[3., 2.], N[2., 1.])
    vl = vertices_list(h)
    # Test Vertices
    @test length(vl) == 4
    @test N[1., 1.] ∈ vl
    @test N[1., 3.] ∈ vl
    @test N[5., 1.] ∈ vl
    @test N[5., 3.] ∈ vl
    # norm
    @test norm(h) == norm(N[5., 3.], Inf)
    # radius
    @test radius(h) == norm(N[2., 1.], Inf)
    # diameter
    @test diameter(h) == norm(N[5., 3.] - N[1., 1.], Inf)

    # alternative constructor
    c = ones(N, 2)
    r = to_N(N, [0.1, 0.2])
    l = to_N(N, [0.9, 0.8])
    h = to_N(N, [1.1, 1.2])
    H1 = Hyperrectangle(c, r)
    H2 = Hyperrectangle(low=l, high=h)
    @test H1.center ≈ H2.center
    @test H1.radius ≈ H2.radius
    @static if VERSION < v"0.7-"
        @test_throws ArgumentError Hyperrectangle(xyz="zyx")
    end

    # Test low and high methods for a hyperrectangle
    H = Hyperrectangle(to_N(N, [-2.1, 5.6, 0.9]), fill(to_N(N, 0.5), 3))
    @test low(H) ≈ to_N(N, [-2.6, 5.1, 0.4])
    @test high(H) ≈ to_N(N, [-1.6, 6.1, 1.4])

    # membership
    H = Hyperrectangle(N[1.0, 1.0], N[2.0, 3.0])
    @test !∈(N[-1.1, 4.1], H)
    @test ∈(N[-1.0, 4.0], H)

    # an_element function
    H = Hyperrectangle(N[1.0, 2.0], N[3.0, 4.0])
    @test an_element(H) ∈ H

    # subset
    H1 = Hyperrectangle(N[1.0, 3.0], N[0.5, 0.5])
    H2 = Hyperrectangle(N[2.0, 2.5], N[0.5, 0.5])
    H3 = Hyperrectangle(N[2.0, 2.0], N[2.0, 3.0])
    B1 = BallInf(N[2.0, 2.5], N(0.5))
    B2 = BallInf(N[2.0, 2.0], N(1.0))
    @test !⊆(H1, H2) && ⊆(H1, H3) && ⊆(H2, H3)
    subset, point = ⊆(H1, H2, true)
    @test !subset && point ∈ H1 && point ∉ H2
    subset, point = ⊆(H2, H1, true)
    @test !subset && point ∈ H2 && point ∉ H1
    subset, point = ⊆(H1, H3, true)
    @test subset
    @test ⊆(H2, B1) && ⊆(B1, H2)
    @test ⊆(B1, B2) && !⊆(B2, B1)

    # intersection & intersection emptiness
    H1 = Hyperrectangle(N[1.0, 1.0], N[2.0, 2.0])
    H2 = Hyperrectangle(N[3.0, 3.0], N[2.0, 2.0])
    B1 = BallInf(N[2.0, 4.0], N(0.5))
    for (X1, X2) in [(H1, H2), (H2, H1)]
        intersection_empty, point = is_intersection_empty(X1, X2, true)
        cap = intersection(X1, X2)
        @test cap isa Hyperrectangle{N} && center(cap) == N[2., 2.] &&
              radius_hyperrectangle(cap) == N[1., 1.]
        @test !is_intersection_empty(X1, X2) &&
              !intersection_empty && point ∈ X1 && point ∈ X2
    end
    cap = intersection(H1, B1)
    @test cap isa EmptySet{N}
    @test is_intersection_empty(H1, B1) && is_intersection_empty(H1, B1, true)[1]

    # linear map (concrete)
    P = linear_map(N[1. 0.; 0. 2.], H1)
    @test P isa VPolygon # in 2D we get a polygon
    P = linear_map(Diagonal(N[1., 2., 3., 4.]),
                   Approximations.overapproximate(H1 * H1))
    @test P isa VPolytope # in 4D we get a VPolytope
end
