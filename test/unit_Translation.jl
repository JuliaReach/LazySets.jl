for N in [Float64, Rational{Int}, Float32]

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
    @test σ(N[1, 0, 0], tr) == N[2, 1, 1]

    # support function
    @test ρ(N[1, 0, 0], tr) == N(2)

    # boundedness
    @test isbounded(tr)
    @test !isbounded(Translation(Universe{N}(3), v))

    if N ∈ [Float64]
        @test !isbounded(Translation(HPolyhedron([HalfSpace(N[1, 1], N(1))]), N[0, 0]))
    end

    # function to get an element
    @test an_element(tr) ∈ tr.X

    # emptiness check
    @test !isempty(tr)

    # set membership
    B = Ball2(zeros(N, 2), N(1))
    @test center(B) ∈ B ⊕ zeros(N, 2)

    # ==================================
    # Type-specific methods
    # ==================================

    # the translation is the origin => constraints remain unchanged
    constraints_list(B ⊕ zeros(N, 3)) == constraints_list(tr.X)

    # one-dimensional case: translation to the left
    Tleft = Translation(HalfSpace(N[1], N(1)), [N(-1)])
    @test constraints_list(Tleft)[1] == HalfSpace(N[1], N(0))

    # one-dimensional case: translation to the right
    Tright = Translation(HalfSpace(N[1], N(1)), [N(1)])
    @test constraints_list(Tright)[1] == HalfSpace(N[1], N(2))

    # translation of a set not represented by a finite number of constraints
    tr = Ball2(zeros(N, 2), N(1)) ⊕ N[1, 0]
    @test ρ(N[1, 0], tr) == N(2)
    @test_throws MethodError constraints_list(tr)
end
