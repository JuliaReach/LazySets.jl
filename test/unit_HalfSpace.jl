for N in [Float64, Rational{Int}, Float32]
    # normal constructor
    normal = ones(N, 3)
    hs = HalfSpace(normal, N(5))

    # dimension
    @test dim(hs) == 3

    # support vector and membership function
    function test_svec(hs, d)
        @test σ(d, hs) ∈ hs
        @test σ(N(2) * d, hs) ∈ hs
        d2 = N[1, 0, 0]
        @test_throws ErrorException σ(d2, hs)
        d2 = zeros(N, 3)
        @test σ(d2, hs) ∈ hs
    end
    # tests 1
    normal = ones(N, 3)
    d = ones(N, 3)
    test_svec(HalfSpace(normal, N(5)), d)
    # tests 2
    normal = zeros(N, 3); normal[3] = N(1)
    d = zeros(N, 3); d[3] = 1
    test_svec(HalfSpace(normal, N(5)), d)

    # check that support vector in opposite direction of normal does not work
    @test_throws ErrorException σ(N[-1], HalfSpace(N[1], N(1)))

    # an_element function and membership function
    @test an_element(hs) ∈ hs

    # constraints list
    @test constraints_list(hs) == [hs]

    # constrained dimensions
    @test constrained_dimensions(HalfSpace(N[1, 0, 1], N(1))) == [1, 3]
    @test constrained_dimensions(HalfSpace(N[0, 1, 0], N(1))) == [2]

    # halfspace_left & halfspace_right
    @test N[1, 2] ∈ halfspace_left(N[1, 1], N[2, 2])
    @test N[2, 1] ∈ halfspace_right(N[1, 1], N[2, 2])

    # intersection emptiness
    b = BallInf(N[3, 3, 3], N(1))
    empty_intersection, v = is_intersection_empty(b, hs, true)
    @test is_intersection_empty(b, hs) && empty_intersection
    b = BallInf(N[1, 1, 1], N(1))
    empty_intersection, v = is_intersection_empty(b, hs, true)
    @test !is_intersection_empty(b, hs) && !empty_intersection && v ∈ hs
    hs1 = HalfSpace(N[1, 0], N(1)) # x <= 1
    hs2 = HalfSpace(N[-1, 0], N(-2)) # x >= 2
    empty_intersection, v = is_intersection_empty(hs1, hs2, true)
    @test is_intersection_empty(hs1, hs2) && empty_intersection && v == N[]
    hs3 = HalfSpace(N[-1, 0], N(-1)) # x >= 1
    empty_intersection, v = is_intersection_empty(hs1, hs3, true)
    @test !is_intersection_empty(hs1, hs3) && !empty_intersection && v == N[1, 0]
    hs4 = HalfSpace(N[-1, 0], N(0)) # x >= 0
    empty_intersection, v = is_intersection_empty(hs1, hs4, true)
    @test !is_intersection_empty(hs1, hs4) && !empty_intersection && v ∈ hs1 && v ∈ hs4
end
