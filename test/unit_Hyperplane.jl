for N in [Float64, Rational{Int}, Float32]
    # random hyperplane
    rand(Hyperplane)

    # normal constructor
    a = ones(N, 3)
    hp = Hyperplane(a, N(5))

    # dimension
    @test dim(hp) == 3

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
    a = ones(N, 3)
    test_svec(Hyperplane(a, N(5)))
    # tests 2
    a = zeros(N, 3); a[3] = N(1)
    test_svec(Hyperplane(a, N(5)))

    # isempty
    @test !isempty(hp)

    # an_element and membership
    @test an_element(hp) ∈ hp

    # constrained dimensions
    @test constrained_dimensions(Hyperplane(N[1, 0, 1], N(1))) == [1, 3]
    @test constrained_dimensions(Hyperplane(N[0, 1, 0], N(1))) == [2]

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
