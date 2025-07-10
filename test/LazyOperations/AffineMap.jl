for N in [Float64, Float32, Rational{Int}]

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

    # linear map of an affine map is automatically simplified to an affine map
    Mam = M * am
    @test Mam isa AffineMap && Mam.M == M * am.M && Mam.v == M * am.v

    # scaling of an affine map is an affine map
    α = N(2)
    αam = α * am
    @test αam isa AffineMap && αam.M == α * am.M && αam.v == α * am.v

    # support vector
    sv = σ(N[1, 0, 0], am)
    @test sv[1] == N(2) &&
          sv[2:3] ∈ linear_map(Diagonal(N[2, 3]), BallInf(zeros(N, 2), N(1)))

    # support function
    @test ρ(N[1, 0, 0], am) == N(2)

    # boundedness
    @test isbounded(am) && isboundedtype(typeof(am))
    am2 = AffineMap(M, Universe{N}(3), v)
    @test !isbounded(am2) && !isboundedtype(typeof(am2))

    # ispolyhedral
    @test ispolyhedral(am)
    if N isa AbstractFloat
        am3 = AffineMap(M, Ball2(zeros(N, 3), N(1)), v)
        @test !ispolyhedral(am3)
    end

    # function to get an element
    @test (an_element(am) - am.v) ∈ (am.M * am.X)
    @test an_element(am) ∈ am.M * am.X ⊕ am.v

    # emptiness check
    @test !isempty(am)
    @test isempty(AffineMap(M, EmptySet{N}(2), v))

    # containment
    L = LineSegment(N[1, 0], N[2, 0])
    b2 = N[1, 0]
    b3 = N[1, 0, 0]
    # well-conditioned square matrix
    M = N[1 2; 3 4]
    @test N[2, 3] ∈ M * L + b2
    @test N[0, 0] ∉ M * L + b2
    # well-conditioned rectangular matrix
    M = N[1 2; 3 4; 5 6]
    @test N[2, 3, 5] ∈ M * L + b3
    @test N[0, 0, 0] ∉ M * L + b3
    # ill-conditioned square matrix
    M = N[-1 -2; 1 2]
    @test N[0, 1] ∈ M * L + b2
    @test N[0, 0] ∉ M * L + b2
    # ill-conditioned rectangular matrix
    M = N[-1 -2; 1 2; 5 6]
    @test N[0, 1, 5] ∈ M * L + b3
    @test N[0, 0, 0] ∉ M * L + b3

    # volume
    B = BallInf(N[0, 0], N(1))
    v = N[-1, 0]
    M = N[1 0; 0 1]
    @test volume(M * B + v) == N(4)
    M = N[1 2; 3 4]
    @test volume(M * B + v) ≈ N(8)
    M = N[-1 -2; -3 -4]
    @test volume(M * B + v) ≈ N(8)
    M = N[0 0; 0 0]
    @test volume(M * B + v) == N(0)
    M = N[0 0;]
    @test_throws DimensionMismatch volume(M * B + N[1])

    # ==================================
    # Type-specific methods
    # ==================================

    # an affine map of the form I*X + b where I is the identity matrix is a pure translation
    #v = N[1, 0, 2]
    #am_tr = AffineMap(I, B, v) # crashes, see #1544
    #@test am_tr isa Translation && am_tr.v == v

    # two-dimensional case
    M = N[1 0; 0 2]
    am = AffineMap(M, B, v)

    # list of vertices check
    vlist = vertices_list(am)
    @test ispermutation(vlist, [N[0, 2], N[-2, 2], N[0, -2], N[-2, -2]])

    # inclusion check
    h = Hyperrectangle(N[-1, 0], N[1, 2])
    @test h ⊆ am && am ⊆ h

    # concretize
    @test concretize(am) == affine_map(M, B, v)
end

for N in [Float64, Float32]
    B = BallInf(zeros(N, 3), N(1))

    # the translation is the origin and the linear map is the identity => constraints remain unchanged
    Id3 = Matrix(one(N) * I, 3, 3)
    @test ispermutation(constraints_list(AffineMap(Id3, B, zeros(N, 3))),
                        constraints_list(B))
end
