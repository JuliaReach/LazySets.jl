for N in [Float64, Rational{Int}, Float32]
    # constructor
    b = BallInf(N[2, 2, 2], N(1))
    s = Singleton(N[0, 0])
    r = Dict(1 => N(4), 3 => N(0))
    r_none = Dict{Int, N}()  # no reset
    r_1 = Dict{Int, N}(1 => N(1))  # reset x1
    r_12 = Dict{Int, N}(1 => N(1), 2 => N(1))  # reset x1 and x2
    rm = ResetMap(b, r)

    # dimension
    @test dim(rm) == 3

    # support vector
    @test σ(N[1, 1, 1], rm) == N[4, 3, 0]
    b1 = Ball1(N[1, 1], N(1))
    @test σ(N[1, 1], ResetMap(b1, Dict(1 => N(3)))) == N[3, 2]
    @test σ(N[1, 1], ResetMap(b1, Dict(2 => N(3)))) == N[2, 3]
    svec = σ(N[1, 1], ResetMap(s, r_1))
    @test svec == N[1, 0] && element(s) == N[0, 0]  # does not modify s

    # boundedness
    @test isbounded(rm)  # bounded set
    @test isbounded(ResetMap(Singleton(N[1, 2]), r_none))  # bounded set without resets
    @test !isbounded(ResetMap(HPolyhedron([HalfSpace(N[1, 1], N(1))]), r_none))  # unbounded set without resets
    @test !isbounded(ResetMap(Universe{N}(2), r_1))  # unbounded set without enough resets
    @test isbounded(ResetMap(Universe{N}(2), r_12))  # unbounded set with enough resets
    if N in [Float64]
        @test !isbounded(ResetMap(HPolyhedron([HalfSpace(N[1, 1], N(1))]), r_1))  # unbounded set without enough resets
        @test isbounded(ResetMap(HPolyhedron([HalfSpace(N[1, 1], N(1))]), r_12))  # unbounded set, but captured by resets
    end

    # an_element function
    an_element(rm)
    elem = an_element(ResetMap(s, r_1))
    @test elem == N[1, 0] && element(s) == N[0, 0]  # does not modify s

    # isempty
    @test !isempty(rm)

    # getter for affine map
    @test get_A(rm) == sparse([2], [2], N[1], 3, 3)
    @test get_b(rm) == sparsevec([1], N[4], 3)

    # constraints_list
    if test_suite_polyhedra
        p = HPolytope(constraints_list(rm))
        @test N[4, 1, 0] ∈ p && N[4, 3, 0] ∈ p && N[2, 2, 2] ∉ p
    end
end
