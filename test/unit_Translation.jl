for N in [Float64, Rational{Int}, Float32]

    # ==================================
    # Constructor and interface methods
    # ==================================
 
    B = BallInf(N[2, 2, 2], N(1))
    v = N[0, 0, 0]
    tr = Translation(B, v)

    # dimension check
    @test_throws AssertionError Translation(B, N[0, 0])
    
    @test dim(tr) == 3

    # support vector
    @test σ(N[1, 1, 1], rm) == N[3, 3, 3]

    # support function
    @test ρ(N[1, 1, 1], tr) == N(9)

    # boundedness
    @test isbounded(tr)
    @test !isbounded(Translation(Universe{N}(3), v))
    
    if N ∈ [Float64]
        @test !isbounded(Translation(HPolyhedron([HalfSpace(N[1, 1], N(1))]), N[0, 0]))
        @test !isbounded(Translation(HPolyhedron([HalfSpace(N[1, 1], N(1))]), N[0, 0]))
    end

    # function to get an element
    @test an_element(tr) ∈ tr.X

    # emptiness check
    @test !isempty(tr)

    # ==================================
    # Type-specific methods
    # ==================================

    # the translation is the origin => constraints remain unchanged
    constraints_list(tr) == constraints_list(tr.X)

    # one-dimensional case: translation to the left
    Tleft = Translation(HalfSpace(N[1], N(1)), [N(-1)])
    @test constraints_list(Tleft)[1] == HalfSpace(N[1], N(0))

    # one-dimensional case: translation to the right
    Tright = Translation(HalfSpace(N[1], N(1)), [N(1)])
    @test constraints_list(Tright)[1] == HalfSpace(N[1], N(2))
end
