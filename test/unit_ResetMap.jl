for N in [Float64, Rational{Int}, Float32]
    # constructor
    b = BallInf(N[2, 2, 2], N(1))
    s = Singleton(N[0, 0])
    r = Dict(1 => N(4), 3 => N(0))
    r_none = Dict{Int, N}()  # no reset
    r_1 = Dict{Int, N}(1 => N(1))  # reset x1
    r_12 = Dict{Int, N}(1 => N(1), 2 => N(1))  # reset x1 and x2
    rm = ResetMap(b, r)
    # rm is the line segment in 3D from [4, 1, 0] to [4, 3, 0]

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
    @test !isbounded(ResetMap(Universe{N}(2), r_1))  # unbounded set without enough resets
    @test isbounded(ResetMap(Universe{N}(2), r_12))  # unbounded set with enough resets

    # an_element function
    an_element(rm)
    elem = an_element(ResetMap(s, r_1))
    @test elem == N[1, 0] && element(s) == N[0, 0]  # does not modify s

    # isempty
    @test !isempty(rm)

    # getter for affine map
    @test matrix(rm) == Diagonal(N[0, 1, 0])
    @test vector(rm) == sparsevec([1], N[4], 3)

    # constraints_list
    if test_suite_polyhedra
        p = HPolytope(constraints_list(rm))
        @test N[4, 1, 0] ∈ p && N[4, 3, 0] ∈ p && N[2, 2, 2] ∉ p
    end
    # constraints_list of a hyperrectangular set
    @test ispermutation(constraints_list(rm),
                        [HalfSpace(N[1, 0, 0], N(4)),
                         HalfSpace(N[-1, 0, 0], N(-4)),
                         HalfSpace(N[0, 1, 0], N(3)),
                         HalfSpace(N[0, -1, 0], N(-1)),
                         HalfSpace(N[0, 0, 1], N(0)),
                         HalfSpace(N[0, 0, -1], N(0))])
end

for N in [Float64]
    b = BallInf(N[2, 2, 2], N(1))
    r = Dict(1 => N(4), 3 => N(0))
    r_none = Dict{Int, N}()  # no reset
    r_1 = Dict{Int, N}(1 => N(1))  # reset x1
    r_12 = Dict{Int, N}(1 => N(1), 2 => N(1))  # reset x1 and x2
    rm = ResetMap(b, r)

    # boundedness
    P = HPolyhedron([HalfSpace(N[1, 1], N(1))])
    @test !isbounded(ResetMap(P, r_none))  # unbounded set without resets
    @test !isbounded(ResetMap(P, r_1))  # unbounded set without enough resets
    @test isbounded(ResetMap(P, r_12))  # unbounded set, but captured by resets
    @test !isbounded(ResetMap(HPolyhedron([HalfSpace(N[1, 1], N(1))]), r_none))  # unbounded set without resets
    @test !isbounded(ResetMap(HPolyhedron([HalfSpace(N[1, 1], N(1))]), r_1))  # unbounded set without enough resets
    @test isbounded(ResetMap(HPolyhedron([HalfSpace(N[1, 1], N(1))]), r_12))  # unbounded set, but captured by resets

    # intersection
    b2 = BallInf(N[4, 2, 0], N(1))
    for cap in [intersection(b2, rm), intersection(rm, b2)]
        @test N[4, 1, 0] ∈ cap && N[4, 3, 0] ∈ cap && N[3, 1, 0] ∉ cap &&
              N[3, 1, 2] ∉ cap
    end
end
