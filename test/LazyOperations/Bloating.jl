for N in [Float64, Float32]
    B = BallInf(zeros(N, 2), N(1))
    E = HPolyhedron([HalfSpace(N[1], N(0)), HalfSpace(N[-1], N(-1))])  # empty
    U = Universe{N}(3)
    ε = N(1 // 10)
    p = N(2)

    # constructors
    X = Bloating(B, ε, p)
    X⁻ = Bloating(B, -ε, p)
    Y = Bloating(E, ε, p)
    Z = Bloating(U, ε, p)
    @test_throws AssertionError Bloating(B, ε, N(9 // 10))

    # dimension
    @test dim(X) == dim(X⁻) == 2 && dim(Y) == 1 && dim(Z) == 3

    # isbounded
    @test isbounded(X) && isbounded(X⁻) && isbounded(Y) && !isbounded(Z) &&
          isboundedtype(typeof(X)) && !isboundedtype(typeof(Z))

    # isempty
    @test !isempty(X) && !isempty(X⁻) && isempty(Y) && !isempty(Z)

    # an_element
    for S in [X, X⁻, Z]
        v = an_element(S)
        @test v isa AbstractVector{N} && length(v) == dim(S)
    end
    @test_throws ErrorException an_element(Y)

    # center
    @test center(X) == center(B)

    # concretize
    B2 = Bloating(B, ε, N(Inf))
    @test concretize(B2) == BallInf(N[0, 0], N(11 // 10))

    # tests for different norms
    for p in N[1, 2, Inf]
        B = BallInf(zeros(N, 2), N(1))
        X = Bloating(B, ε, p)
        X⁻ = Bloating(B, -ε, p)
        E = HPolyhedron([HalfSpace(N[1], N(0)), HalfSpace(N[-1], N(-1))])
        Y = Bloating(E, ε, p)  # empty set
        Bp = Ballp(p, zeros(N, 2), ε)

        # support vector and support function
        d = N[1, 0]
        @test σ(d, X) == σ(d, B + Bp)
        @test ρ(d, X) == ρ(d, B + Bp) == N(1 + ε)
        @test ρ(d, X⁻) == N(1 - ε)
        if p != N(1)
            @test ρ(d, X⁻) == dot(d, σ(d, X⁻))
        end
        @test_throws ArgumentError σ(N[1], Y)
        @test_throws ArgumentError ρ(N[1], Y)
        d = N[1, 99 // 100]
        @test σ(d, X) == σ(d, B + Bp)
        @test ρ(d, X) == ρ(d, B + Bp)

        # subset
        if p == N(Inf)
            @test X⁻ ⊆ B ⊆ X
        elseif p == N(1)
            @test B ⊆ X
        end

        # bloating of negative bloating gives the identity
        if p != N(1)
            # need to overapproximate because constraints_list is missing
            B2 = box_approximation(Bloating(Bloating(B, ε, p), -ε, p))
            @test isequivalent(B, B2)
        end

        # tests for infinity norm
        if p == Inf
            B_bloated = BallInf(zeros(N, 2), N(11 / 10))
            @test constraints_list(X) == constraints_list(B_bloated)
        end
    end
end
