for N in [Float64, Rational{Int}, Float32]
    # random hyperplane
    rand(Hyperplane)

    # normal constructor
    a = ones(N, 3)
    b = N(5)
    hp = Hyperplane(a, b)

    # corner case: zero normal vector
    @test_throws AssertionError Hyperplane(N[0, 0], N(1))

    # dimension
    @test dim(hp) == 3

    # support function
    hp2 = Hyperplane(N[1, 0], N(1))
    @test ρ(N[2, 0], hp2) == N(2)
    @test ρ(N[-2, 0], hp2) == N(-2)
    @test ρ(N[1, 1], hp2) == N(Inf)

    # support vector and membership
    function test_svec(hp)
        d1 = copy(hp.a)
        @test σ(d1, hp) ∈ hp
        @test σ(N(2) * d1, hp) ∈ hp
        d2 = N[1, 0, 0]
        @test_throws ErrorException σ(d2, hp)
        d3 = N[1, 1, 2]
        @test_throws ErrorException σ(d3, hp)
        d4 = zeros(N, 3)
        @test σ(d4, hp) ∈ hp
    end
    # tests 1
    test_svec(hp)
    # tests 2
    test_svec(Hyperplane(N[0, 0, 1], b))

    # support vector in opposite direction
    hp2 = Hyperplane(N[1], N(1))
    @test σ(N[-1], hp2) ∈ hp2
    # support vector in other directions throws an error (but see #750)
    @test_throws ErrorException σ(N[1, 1], Hyperplane(N[1, 0], N(1)))

    # boundedness
    @test !isbounded(hp)

    # universality
    @test !isuniversal(hp)
    res, w = isuniversal(hp, true)
    @test !res && w ∉ hp

    # isempty
    @test !isempty(hp)

    # an_element and membership
    @test an_element(hp) ∈ hp

    # constrained dimensions
    @test constrained_dimensions(Hyperplane(N[1, 0, 1], N(1))) == [1, 3]
    @test constrained_dimensions(Hyperplane(N[0, 1, 0], N(1))) == [2]

    # constraints_list
    @test ispermutation(constraints_list(hp),
                        [HalfSpace(a, b), HalfSpace(-a, -b)])

    # test concrete linear map of a hyperplane
    H = Hyperplane(N[1, -1], N(0)) # x = y
    M = N[1 0; 0 0] # non-invertible matrix
    if N == Float32 || N == Float64
        @test linear_map(M, H) isa Hyperplane{Float64}
    elseif N == Rational{Int}
        @test linear_map(M, H) isa Hyperplane{Rational{BigInt}}
    end
    @test_throws ArgumentError linear_map(M, H, algorithm="inv")
    M = N[2 2; 0 1] # invertible matrix
    @test linear_map(M, H) == Hyperplane(N[0.5, -2.0], N(0.0))

    # translation
    @test translate(hp, N[1, 2, 3]) == Hyperplane(ones(N, 3), N(11))

    # intersection emptiness
    u = Universe{N}(3)
    res, v = is_intersection_empty(u, hp, true)
    @test !is_intersection_empty(u, hp) && !res && v ∈ hp && v ∈ u
    res, v = is_intersection_empty(hp, u, true)
    @test !is_intersection_empty(hp, u) && !res && v ∈ hp && v ∈ u
end

# Polyhedra tests that only work with Float64
for N in [Float64]
    hp = Hyperplane(ones(N, 3), N(5))

    # intersection emptiness
    b = BallInf(zeros(N, 3), N(1))
    empty_intersection, v = is_intersection_empty(b, hp, true)
    @test empty_intersection && is_intersection_empty(b, hp) && v == N[]
    b = BallInf(N[3, 3, 3], N(1))
    empty_intersection, v = is_intersection_empty(b, hp, true)
    @test empty_intersection && is_intersection_empty(b, hp) && v == N[]
    b = BallInf(N[2, 2, 2], N(1))
    empty_intersection, v = is_intersection_empty(b, hp, true)
    @test !empty_intersection && !is_intersection_empty(b, hp) && v ∈ hp && v ∈ b
end
