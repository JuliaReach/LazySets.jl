for N in [Float64, Rational{Int}, Float32]
    # random empty set
    rand(EmptySet)

    E = EmptySet{N}(2)
    B = BallInf(ones(N, 2), N(1))

    # testing that the empty set is an absorbing element for the cartesian product
    @test B * E isa EmptySet && E * B isa EmptySet
    # testing the mathematical alias ×
    @test B × E isa EmptySet && E × B isa EmptySet

    cpa = CartesianProductArray([B, N(2) * B, N(3) * B])
    @test cpa * E isa EmptySet && E * cpa isa EmptySet
    @test cpa × E isa EmptySet && E × cpa isa EmptySet

    # testing cp of empty set with itself
    @test E * E == E

    # testing that the empty set is an absorbing element for the Minkowski sum
    @test B + E isa EmptySet && E + B isa EmptySet
    # testing the mathematical alias ⊕
    @test B ⊕ E isa EmptySet && E ⊕ B isa EmptySet

    msa = MinkowskiSumArray([B, N(2) * B, N(3) * B])
    @test msa + E isa EmptySet && E + msa isa EmptySet
    @test msa ⊕ E isa EmptySet && E ⊕ msa isa EmptySet

    # testing M-sum of empty set with itself
    @test E + E == E

    # testing that the emptyset is neutral for the convex hull
    @test CH(B, E) == B
    @test CH(E, B) == B

    # test convex hull of empty set with itself
    @test CH(E, E) == E

    # dim
    @test dim(E) == 2

    # support vector
    @test_throws ErrorException σ(N[0], E)

    # boundedness
    @test isbounded(E)

    # isuniversal
    res, w = isuniversal(E, true)
    @test !isuniversal(E) && !res && w ∉ E

    # membership
    @test N[0] ∉ E
    @test N[0, 0] ∉ E

    # subset
    @test E ⊆ B && ⊆(E, B, true)[1]
    subset, point = ⊆(B, E, true)
    @test B ⊈ E && !subset && point ∈ B
    @test E ⊆ E && ⊆(E, E, true)[1]

    # emptiness check
    @test isempty(E)

    # an_element/norm/radius/diameter functions
    @test_throws ErrorException an_element(E)
    @test_throws ErrorException norm(E)
    @test_throws ErrorException radius(E)
    @test_throws ErrorException diameter(E)

    # linear map of an empty set
    linear_map(ones(N, 2, 2), E) == E

    # translation
    @test translate(E, N[1, 2]) == E
end

# default Float64 constructor
@test EmptySet(2) == ∅(2) == EmptySet{Float64}(2)
