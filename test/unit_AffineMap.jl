for N in [Float64, Rational{Int}, Float32]

    # ==================================
    # Constructor and interface methods
    # ==================================
 
    B = BallInf(zeros(N, 3), N(1))
    v = N[1, 0, 0] # translation along dimension 1
    M = Diagonal([1.0, 2.0, 3.0])
    am = AffineMap(M, v, B)

    # dimension check
    @test dim(am) == 3

    # dimension assertion
    @test_throws AssertionError AffineMap(M, N[0, 0], B)
    @test_throws AssertionError AffineMap(M, v, B × B)

    # linear map of an affine map is automatically simplified to a linear map
    Mam = M * am
    @test Mam isa AffineMap && Mam.M == M * am.M && Mam.b == M * am.b

    # scaling of an affine map is an affine map
    α = 2.0
    αam = α * am
    @test αam isa AffineMap && αam.M == α * am.M && αam.b == α * am.b

    # support vector
    @test σ(N[1, 0, 0], am) == N[2, 2, 3]

    # support function
    @test ρ(N[1, 0, 0], am) == N(2)

    # boundedness
    @test isbounded(am)
    @test !isbounded(AffineMap(M, v, Universe{N}(3)))

    # function to get an element
    @test (an_element(am) - am.b) ∈ (am.M * am.X) # an_element(am) ∈ am.M * am.X ⊕ am.b : requires #1358

    # emptiness check
    @test !isempty(am)

    # ==================================
    # Type-specific methods
    # ==================================

    # the translation is the origin and the linear map is the identity => constraints remain unchanged
    Id3 = Matrix(one(N) * I, 3, 3)
    @test constraints_list(AffineMap(Id3, zeros(N, 3), B)) == constraints_list(B)

    # two-dimensional case: inclusion check
    B2 = BallInf(zeros(N, 2), N(1))
    M = N[1.0 0.0; 0.0 2.0]
    tl = N[-1.0, 0.0]
    am = AffineMap(M, tl, B2)

    # list of vertice check
    vlist = vertices_list(am)
    @test ispermutation(vlist, [N[0, 2], N[-2, 2], N[0, -2], N[-2, -2]])

    # inclusion check
    h = Hyperrectangle(N[-1, 0], N[1, 2])
    @test h ⊆ am && am ⊆ h
end
