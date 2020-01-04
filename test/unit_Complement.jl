for N in [Float64, Rational{Int}, Float32]
    B1 = BallInf(N[0, 0], N(1))
    B2 = BallInf(N[4, -4], N(1))
    B3 = BallInf(N[1, -1], N(1))
    C = Complement(B1)

    # double-complement is the identity
    @test Complement(C) == B1

    # dimension
    @test dim(C) == 2

    # membership
    @test N[1, 1] ∉ C && N[2, 2] ∈ C

    # emptiness
    @test !isempty(C) && isempty(Complement(Universe(2)))

    # inclusion
    subset, point = ⊆(B2, C, true)
    @test B2 ⊆ C && subset && point == N[]
    for X in [B1, B3]
        subset, point = ⊆(X, C, true)
        @test X ⊈ C && !subset && point ∈ X && point ∉ C
    end

    # isdisjoint
    for X in [B2, B3]
        res, w = isdisjoint(X, C, true)
        @test !isdisjoint(X, C) && !res && w ∈ X && w ∈ C
        res, w = isdisjoint(C, X, true)
        @test !isdisjoint(C, X) && !res && w ∈ X && w ∈ C
    end
    res, w = isdisjoint(B1, C, true)
    @test isdisjoint(B1, C) && res && w == N[]
    res, w = isdisjoint(C, B1, true)
    @test isdisjoint(C, B1) && res && w == N[]

    # test convexity from the type
    @test isconvextype(typeof(Complement(Universe{N}(2))))
    @test isconvextype(typeof(Complement(EmptySet{N}(2))))
    @test isconvextype(typeof(Complement(HalfSpace(N[1], N(0)))))
end
