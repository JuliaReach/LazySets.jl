for N in [Float64, Float32]
    B = BallInf(zeros(N, 2), N(1))
    E = HPolyhedron([HalfSpace(N[1], N(0)), HalfSpace(N[-1], N(-1))])  # empty
    U = Universe{N}(3)
    ε = N(1//10)
    p = N(2)

    # constructors
    X = Bloating(B, ε, p)
    Y = Bloating(E, ε, p)
    Z = Bloating(U, ε, p)
    @test_throws AssertionError Bloating(B, N(0), p)
    @test_throws AssertionError Bloating(B, ε, N(9//10))

    # dimension
    @test dim(X) == 2 && dim(Y) == 1 && dim(Z) == 3

    # isbounded
    @test isbounded(X) && !isbounded(Z)
    # @test isbounded(Y)  # currently crashes (see #1779)

    # isempty
    @test !isempty(X) && isempty(Y) && !isempty(Z)

    # an_element
    v = an_element(X)
    @test v isa AbstractVector{N} && length(v) == 2
    @test_throws ErrorException an_element(Y)
    v = an_element(Z)
    @test v isa AbstractVector{N} && length(v) == 3

    # tests for different norms
    for p in N[1, 2, Inf]
        B = BallInf(zeros(N, 2), N(1))
        X = Bloating(B, ε, p)
        E = HPolyhedron([HalfSpace(N[1], N(0)), HalfSpace(N[-1], N(-1))])
        Y = Bloating(E, ε, p)  # empty set
        Bp = Ballp(p, zeros(N, 2), ε)

        # support vector and support function
        d = N[1, 0]
        @test σ(d, X) == σ(d, B + Bp)
        @test ρ(d, X) == ρ(d, B + Bp) == N(1 + ε)
        @test_throws ErrorException σ(N[1], Y)
        @test_throws ErrorException ρ(N[1], Y)
        d = N[1, 99//100]
        @test σ(d, X) == σ(d, B + Bp)
        @test ρ(d, X) == ρ(d, B + Bp)
    end
end
