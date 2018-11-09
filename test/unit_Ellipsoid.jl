for N in [Float64, Float32]
    # random ellipsoid
    rand(Ellipsoid)

    # 1D ellipsoid
    E = Ellipsoid(N[0], Diagonal(N[1]))
    # Test Dimension
    @test dim(E) == 1
    # Test Support Vector
    d = N[1]
    @test σ(d, E) == N[1]
    d = N[-1]
    @test σ(d, E) == N[-1]
    # test constructor
    E = Ellipsoid(Diagonal(N[1]))
    @test E.center == N[0]

    # 2D Ellipsoid
    E = Ellipsoid(N[0, 0], Diagonal(N[1, 1]))
    # Test Dimension
    @test dim(E) == 2
    # Test Support Vector
    d = N[1, 0]
    @test σ(d, E) == N[1, 0]
    d = N[-1, 0]
    @test σ(d, E) == N[-1, 0]
    d = N[0, 1]
    @test σ(d, E) == N[0, 1]
    d = N[0, -1]
    @test σ(d, E) == N[0, -1]
    d = N[0, 0]
    @test σ(d, E) ∈ E

    # 2D Ellipsoid not 0-centered
    E = Ellipsoid(N[1, 2], Diagonal(N[1, 1]))
    # Test Dimension
    @test dim(E) == 2
    # Test Support Vector
    d = N[1, 0]
    @test σ(d, E) == N[2, 2]
    d = N[-1, 0]
    @test σ(d, E) == N[0, 2]
    d = N[0, 1]
    @test σ(d, E) == N[1, 3]
    d = N[0, -1]
    @test σ(d, E) == N[1, 1]

    # another shape matrix
    E = Ellipsoid(N[1, 2], Diagonal(N[0.5, 2]))
    # Test Support Vector
    d = N[1, 0]
    @test σ(d, E) ≈ N[1+sqrt(0.5), 2]
    d = N[-1, 0]
    @test σ(d, E) ≈ N[1-sqrt(0.5), 2]
    d = N[0, 1]
    @test σ(d, E) ≈ N[1, 2 + sqrt(2)]
    d = N[0, -1]
    @test σ(d, E) ≈ N[1, 2 - sqrt(2)]

    # isempty
    @test !isempty(E)

    # an_element and set membership functions
    E = Ellipsoid(N[1, 2], Matrix{N}(2I, 2, 2))
    @test an_element(E) ∈ E
end
