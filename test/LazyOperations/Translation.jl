for N in @tN([Float64, Float32, Rational{Int}])
    # ==================================
    # Constructor and interface methods
    # ==================================

    B = BallInf(zeros(N, 3), N(1))
    v = N[1, 0, 0] # translation along dimension 1
    tr = Translation(B, v)

    # test alternative constructors
    @test tr == B + v == B ⊕ v
    # check that translation from the left works as well
    @test tr == v + B == v ⊕ B

    # dimension check
    @test_throws AssertionError Translation(B, N[0, 0])

    @test dim(tr) == 3

    # support vector
    sv = σ(N[1, 0, 0], tr)
    @test sv[1] == N(2) && sv[2:3] ∈ BallInf(zeros(N, 2), N(1))

    # support function
    @test ρ(N[1, 0, 0], tr) == N(2)

    # boundedness
    @test isbounded(tr) && isboundedtype(typeof(tr))
    tr2 = Translation(Universe{N}(3), v)
    @test !isbounded(tr2) && !isboundedtype(typeof(tr2))

    if N ∈ [Float64]
        @test !isbounded(Translation(HPolyhedron([HalfSpace(N[1, 1], N(1))]), N[0, 0]))
    end

    # ispolyhedral
    @test ispolyhedral(tr)
    if N isa AbstractFloat
        tr2 = Translation(Ball2(zeros(N, 3), N(1)), v)
        @test !ispolyhedral(tr2)
    end

    # function to get an element
    @test an_element(tr) ∈ tr.X

    # emptiness check
    @test !isempty(tr)

    # set membership
    B = Ball1(zeros(N, 2), N(1))
    @test center(B) ∈ B ⊕ zeros(N, 2)

    # concrete linear map
    B = BallInf(ones(N, 2), N(1))
    v = N[1, 0]
    tr = B ⊕ v
    clm = linear_map(N[1 2; 0 1], tr)
    @test ispermutation(vertices_list(clm), [N[5, 2], [1, 0], [3, 0], [7, 2]])

    # concretize
    @test concretize(tr) == translate(B, v)

    # center
    @test center(tr) == center(B) + v

    # isuniversal
    @test !isuniversal(tr)
    @test isuniversal(Translation(Universe{N}(2), v))

    # ==================================
    # Type-specific methods
    # ==================================

    # the translation is the origin => constraints remain unchanged
    constraints_list(B ⊕ zeros(N, 2)) == constraints_list(tr.X)

    # one-dimensional case: translation to the left
    Tleft = Translation(HalfSpace(N[1], N(1)), [N(-1)])
    @test constraints_list(Tleft)[1] == HalfSpace(N[1], N(0))

    # one-dimensional case: translation to the right
    Tright = Translation(HalfSpace(N[1], N(1)), [N(1)])
    @test constraints_list(Tright)[1] == HalfSpace(N[1], N(2))

    # the translation of a lazy linear map returns an affine map
    M = N[1 0; 0 2]
    B = BallInf(zeros(N, 2), N(1))
    v = N[1, 0]
    tr = M * B ⊕ v
    @test tr isa AffineMap && tr.M == M && tr.X == B && tr.v == v
end

for N in @tN([Float64, Float32])
    # translation of a set not represented by a finite number of constraints
    tr = Ball2(zeros(N, 2), N(1)) ⊕ N[1, 0]
    @test ρ(N[1, 0], tr) == N(2)
    @test_throws MethodError constraints_list(tr)
end
