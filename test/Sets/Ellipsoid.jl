for N in [Float64, Float32]
    # random ellipsoid
    rand(Ellipsoid)

    # 1D ellipsoid
    E = Ellipsoid(N[0], Diagonal(N[1]))
    # Test Dimension
    @test dim(E) == 1
    # support function
    @test ρ(N[1], E) == ρ(N[-1], E) == N(1)
    # Test Support Vector
    d = N[1]
    @test σ(d, E) == N[1]
    d = N[-1]
    @test σ(d, E) == N[-1]
    # test constructor
    E = Ellipsoid(Diagonal(N[1]))
    @test E.center == center(E) == N[0]
    @test shape_matrix(E) == Diagonal(N[1])

    # 2D Ellipsoid
    E = Ellipsoid(N[0, 0], Diagonal(N[1, 1]))
    # Test Dimension
    @test dim(E) == 2
    # support function
    @test ρ(N[1, 0], E) == ρ(N[-1, 0], E) == N(1)
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

    # boundedness
    @test isbounded(E)

    # isempty
    @test !isempty(E)

    # isuniversal
    answer, w = isuniversal(E, true)
    @test !isuniversal(E) && !answer && w ∉ E

    # an_element and set membership functions
    M = Matrix{N}(2I, 2, 2)
    E = Ellipsoid(N[1, 2], M)
    @test an_element(E) ∈ E

    # translation
    @test translate(E, N[1, 2]) == Ellipsoid(N[2, 4], M)

    # check_posdef optional constructor flag
    @test_throws ArgumentError Ellipsoid(N[1 1; 1 1])
    E = Ellipsoid(N[1 1; 1 1], check_posdef=false)
    @test E.center == zeros(2) && E.shape_matrix == N[1 1; 1 1]

end
