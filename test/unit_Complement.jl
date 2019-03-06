for N in [Float64, Rational{Int}, Float32]
    B1 = BallInf(zeros(N, 2), N(1))
    B2 = BallInf(N[4, -4], N(1))
    B3 = BallInf(N[1, -1], N(1))
    B4 = BallInf(N[0, 0], N(1))
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
    for X in [B1, B3, B4]
        subset, point = ⊆(X, C, true)
        @test !(X ⊆ C) && !subset && point ∈ X && point ∉ C
    end
end
