for N in [Float64, Rational{Int}, Float32]
    E = EmptySet{N}()
    B = BallInf(ones(N, 2), N(1.))

    # testing that the empty set is an absorbing element for the cartesian product
    @test B * E isa EmptySet && E * B isa EmptySet
    # testing the mathematical alias ×
    @test B × E isa EmptySet && E × B isa EmptySet

    cpa = CartesianProductArray([B, N(2.) * B, N(3.) * B])
    @test cpa * E isa EmptySet && E * cpa isa EmptySet
    @test cpa × E isa EmptySet && E × cpa isa EmptySet

    # testing cp of empty set with itself
    @test E * E == E

    # testing that the empty set is an absorbing element for the Minkowski sum
    @test B + E isa EmptySet && E + B isa EmptySet
    # testing the mathematical alias ⊕
    @test B ⊕ E isa EmptySet && E ⊕ B isa EmptySet

    msa = MinkowskiSumArray([B, N(2.) * B, N(3.) * B])
    @test msa + E isa EmptySet && E + msa isa EmptySet
    @test msa ⊕ E isa EmptySet && E ⊕ msa isa EmptySet

    # testing M-sum of empty set with itself
    @test E + E == E

    # testing that the emptyset is neutral for the convex hull
    @test CH(B, E) == B
    @test CH(E, B) == B

    # test convex hull of empty set with itself
    @test CH(E, E) == E

    # an_element function
    @test_throws ErrorException an_element(E)
end

# default Float64 constructor
@test ∅ == EmptySet{Float64}()
