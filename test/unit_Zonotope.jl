for N in [Float64, Rational{Int}, Float32]
    # random zonotope
    rand(Zonotope)

    # 1D Zonotope
    z = Zonotope(N[0], Matrix{N}(I, 1, 1))
    # Test Dimension
    @test dim(z) == 1
    # Test Support Vector
    d = N[1]
    @test σ(d, z) == N[1]
    d = N[-1]
    @test σ(d, z) == N[-1]

    # 2D Zonotope
    z = Zonotope(N[0, 0], Matrix{N}(I, 2, 2))
    # Test Dimension
    @test dim(z) == 2
    # Test Support Vector
    d = N[1, 0]
    @test σ(d, z) == N[1, 1] || N[1, -1]
    d = N[-1, 0]
    @test σ(d, z) == N[-1, 1] || N[-1, -1]
    d = N[0, 1]
    @test σ(d, z) == N[1, 1] || N[-1, 1]
    d = N[0, -1]
    @test σ(d, z) == N[1, -1] || N[-1, -1]

    # 2D Zonotope not 0-centered
    z = Zonotope(N[1, 2], Matrix{N}(I, 2, 2))
    # Test Dimension
    @test dim(z) == 2
    # Test Support Vector
    d = N[1, 0]
    @test σ(d, z) == N[2, 3]
    d = N[-1, 0]
    @test σ(d, z) == N[0, 3]
    d = N[0, 1]
    @test σ(d, z) == N[2, 3]
    d = N[0, -1]
    @test σ(d, z) == N[2, 1]

    # isempty
    @test !isempty(z)

    # an_element function
    @test an_element(z) ∈ z

    # concrete operations
    Z1 = Zonotope(N[1, 1], N[1 1; -1 1])
    Z2 = Zonotope(N[-1, 1], Matrix{N}(I, 2, 2))
    A = N[0.5 1; 1 0.5]

    # concrete Minkowski sum
    Z3 = minkowski_sum(Z1, Z2)
    @test Z3.center == N[0, 2]
    @test Z3.generators == N[1 1 1 0; -1 1 0 1]

    # concrete linear map and scale
    Z4 = linear_map(A, Z3)
    @test Z4.center == N[2, 1]
    @test Z4.generators == N[-0.5 1.5 0.5 1; 0.5 1.5 1 0.5]
    Z5 = scale(0.5, Z3)
    @test Z5.center == N[0, 1]
    @test Z5.generators == N[0.5 0.5 0.5 0; -0.5 0.5 0 0.5]

    # intersection with a hyperplane
    H1 = Hyperplane(N[1, 1], N(3))
    intersection_empty, point = is_intersection_empty(Z1, H1, true)
    @test !is_intersection_empty(Z1, H1) && !intersection_empty &&
          point ∈ Z1 && point ∈ H1
    H2 = Hyperplane(N[1, 1], N(-11))
    @test is_intersection_empty(Z1, H2) && is_intersection_empty(Z1, H2, true)[1]
    @test !is_intersection_empty(H1, Z1)
    @test is_intersection_empty(H2, Z1) && is_intersection_empty(H2, Z1, true)[1]

    # test number of generators
    Z = Zonotope(N[2, 1], N[-0.5 1.5 0.5 1; 0.5 1.5 1 0.5])
    @test ngens(Z) == 4
    # test order reduction
    Zred1 = reduce_order(Z, 1)
    @test ngens(Zred1) == 2
    @test order(Zred1) == 1
    Zred2 = reduce_order(Z, 2)
    @test ngens(Zred2) == 4
    @test order(Zred2) == 2
    Z = Zonotope(N[2, 1], N[-0.5 1.5 0.5 1 0 1; 0.5 1.5 1 0.5 1 0])
    Zred3 = reduce_order(Z, 2)
    @test ngens(Zred3) == 4
    @test order(Zred3) == 2

    # test conversion from hyperrectangular sets
    Z = convert(Zonotope, Hyperrectangle(N[2, 3], N[4, 5]))
    @test Z.center == N[2, 3] && diag(Z.generators) == N[4, 5]
    convert(Zonotope, BallInf(N[5, 3], N(2)))
end
