for N in [Float64, Rational{Int}, Float32]

    # ==================================
    # Constructor and interface methods
    # ==================================
 
    B = BallInf(zeros(N, 3), N(1))
    v = N[1, 0, 0] # translation along dimension 1
    M = Diagonal(N[1, 2, 3])
    am = AffineMap(M, B, v)

    # dimension check
    @test dim(am) == 3

    # dimension assertion
    @test_throws AssertionError AffineMap(M, B, N[0, 0])
    @test_throws AssertionError AffineMap(M, B × B, v)

    # linear map of an affine map is automatically simplified to a linear map
    Mam = M * am
    @test Mam isa AffineMap && Mam.M == M * am.M && Mam.v == M * am.v

    # scaling of an affine map is an affine map
    α = N(2)
    αam = α * am
    @test αam isa AffineMap && αam.M == α * am.M && αam.v == α * am.v

    # support vector
    @test σ(N[1, 0, 0], am) == N[2, 2, 3]

    # support function
    @test ρ(N[1, 0, 0], am) == N(2)

    # boundedness
    @test isbounded(am)
    @test !isbounded(AffineMap(M, Universe{N}(3), v))

    # function to get an element
    @test (an_element(am) - am.v) ∈ (am.M * am.X)
    @test an_element(am) ∈ am.M * am.X ⊕ am.v

    # emptiness check
    @test !isempty(am)

    # ==================================
    # Type-specific methods
    # ==================================

    # the translation is the origin and the linear map is the identity => constraints remain unchanged
    Id3 = Matrix(one(N) * I, 3, 3)
    @test constraints_list(AffineMap(Id3, B, zeros(N, 3))) == constraints_list(B)

    # two-dimensional case
    B2 = BallInf(zeros(N, 2), N(1))
    M = N[1 0; 0 2]
    v = N[-1, 0]
    am = AffineMap(M, B2, v)

    # list of vertices check
    vlist = vertices_list(am)
    @test ispermutation(vlist, [N[0, 2], N[-2, 2], N[0, -2], N[-2, -2]])

    # inclusion check
    h = Hyperrectangle(N[-1, 0], N[1, 2])
    @test h ⊆ am && am ⊆ h
end
