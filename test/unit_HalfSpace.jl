for N in [Float64, Rational{Int}, Float32]
    # random half-space
    rand(HalfSpace)

    # normal constructor
    hs = HalfSpace(ones(N, 3), N(5))

    # corner case: zero normal vector
    @test_throws AssertionError HalfSpace(N[0, 0], N(1))

    # dimension
    @test dim(hs) == 3

    # support function
    hs2 = HalfSpace(N[1, 0], N(1))
    @test ρ(N[2, 0], hs2) == N(2)
    @test ρ(N[-2, 0], hs2) == N(Inf)
    @test ρ(N[1, 1], hs2) == N(Inf)

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

    # support vector in other directions throws an error (but see #750)
    # opposite direction
    @test_throws ErrorException σ(N[-1], HalfSpace(N[1], N(1)))
    # any other direction
    @test_throws ErrorException σ(N[1, 1], HalfSpace(N[1, 0], N(1)))

    # boundedness
    @test !isbounded(hs)

    # universality
    @test !isuniversal(hs)
    res, w = isuniversal(hs, true)
    @test !res && w ∉ hs

    # isempty
    @test !isempty(hs)

    # an_element function and membership function
    @test an_element(hs) ∈ hs

    # constraints list
    @test constraints_list(hs) == [hs]

    # constraints list from matrix-vector representation
    A = N[2 0; 1 3]
    b = N[-1, 1]
    @test constraints_list(A, b) ==
          [HalfSpace(N[2, 0], N(-1)), HalfSpace(N[1, 3], N(1))]

    # constrained dimensions
    @test constrained_dimensions(HalfSpace(N[1, 0, 1], N(1))) == [1, 3]
    @test constrained_dimensions(HalfSpace(N[0, 1, 0], N(1))) == [2]
    # sparse vector
    @test constrained_dimensions(HalfSpace(sparsevec([2], N[1], 3), N(1))) == [2]

    # halfspace_left & halfspace_right
    @test N[1, 2] ∈ halfspace_left(N[1, 1], N[2, 2])
    @test N[2, 1] ∈ halfspace_right(N[1, 1], N[2, 2])

    # translation
    @test translate(hs, N[1, 2, 3]) == HalfSpace(ones(N, 3), N(11))

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

    # check for tighter constraint
    c1 = HalfSpace([1., 0.], 1.)
    c2 = HalfSpace([2., 0.], 2.)
    c3 = HalfSpace([2., 0.], 1.)
    @test LazySets.is_tighter_same_dir_2D(c1, c2) &&
          LazySets.is_tighter_same_dir_2D(c3, c2)
    @test !LazySets.is_tighter_same_dir_2D(c1, c2, strict=true) &&
          LazySets.is_tighter_same_dir_2D(c3, c2, strict=true)

    # test concrete linear map of a half-space
    H = HalfSpace(N[1, -1], N(0)) # x <= y
    M = N[1 0; 0 0] # non-invertible matrix
    @test_throws ArgumentError linear_map(M, H, algorithm="vrep")
    M = N[2 2; 0 1] # invertible matrix
    @test linear_map(M, H) == HalfSpace(N[0.5, -2.0], N(0.0))
end
