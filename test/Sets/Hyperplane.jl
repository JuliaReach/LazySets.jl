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
    @test isbounded(Hyperplane(ones(N, 1), N(1)))

    # is_polyhedral
    @test is_polyhedral(hp)

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
    # projection is y = 0
    if test_suite_polyhedra
        lm = linear_map(M, H)
        if N == Float32 || N == Float64
            @test lm isa Hyperplane{Float64}
            @test lm.a ≈ N[0, -1] && lm.b ≈ N(0)

            # returned set is universal
            @test linear_map(N[1 1], H) == Universe{N}(1)
        elseif N == Rational{Int}
            @test lm isa Hyperplane{Rational{BigInt}, Vector{Rational{BigInt}}}
            @test lm.a == N[0//1, -1//1] && lm.b == N(0//1)
        end
    end

    # projection
    H = Hyperplane(N[1, -1], N(0))  # x = y
    @test project(H, [1]) == project(H, [2]) == Universe{N}(1)
    @test project(H, [1, 2]) == H

    @test_throws ArgumentError linear_map(M, H, algorithm="inv")
    M = N[2 2; 0 1] # invertible matrix

    if test_suite_polyhedra
        @test linear_map(M, H) == Hyperplane(N[0.5, -2.0], N(0.0))
    end

    # translation
    @test translate(hp, N[1, 2, 3]) == Hyperplane(ones(N, 3), N(11))

    # intersection emptiness
    hp2 = Hyperplane(hp.a, hp.b)
    res, v = is_intersection_empty(hp, hp2, true)
    @test !is_intersection_empty(hp, hp2) && !res && v ∈ hp && v ∈ hp2
    hp2 = Hyperplane(hp.a, hp.b + N(1))
    res, v = is_intersection_empty(hp, hp2, true)
    @test is_intersection_empty(hp, hp2) && res && v == N[]
    u = Universe{N}(3)
    res, v = is_intersection_empty(u, hp, true)
    @test !is_intersection_empty(u, hp) && !res && v ∈ hp && v ∈ u
    res, v = is_intersection_empty(hp, u, true)
    @test !is_intersection_empty(hp, u) && !res && v ∈ hp && v ∈ u

    # conversion from polyhedron
    for (P, eq) in [
            (HPolyhedron([HalfSpace(N[-1], N(0)), HalfSpace(N[1//2], N(0))]), true),
            (HPolyhedron([HalfSpace(N[1, 3], N(1)), HalfSpace(N[-2, -6], N(-2))]), true),
            (HPolyhedron([HalfSpace(N[1, 3], N(1)), HalfSpace(N[1, 3], N(1)),
                          HalfSpace(N[-2, -6], N(-2))]), true),
            (HPolyhedron([HalfSpace(N[1, 3], N(1))]), false),
            (HPolyhedron([HalfSpace(N[1, 3], N(1)), HalfSpace(N[1, 4], N(1)),
                          HalfSpace(N[-2, -6], N(-2))]), false),
            (HPolyhedron([HalfSpace(N[1, 3], N(1)), HalfSpace(N[-2, 6], N(-2))]), false),
            (HPolyhedron([HalfSpace(N[1, 3], N(1)), HalfSpace(N[-2, -6], N(2))]), false)
           ]
        @test is_hyperplanar(P) == eq
        if eq
            H = convert(Hyperplane, P)
            @test is_hyperplanar(H)
            @test H isa Hyperplane{N} && isequivalent(P, H)
        else
            @test_throws ArgumentError convert(Hyperplane, P)
        end
    end

    # reflection
    H = Hyperplane(N[1, 1], N(1))
    p = N[0, 0]
    @test reflect(p, H) == N[1, 1]

    # projecting a point onto a line
    H = Hyperplane(N[1, -1], N(0))  # x = y
    @test project(N[1, 0], H) ≈ N[1//2, 1//2]
end

# tests that only work with Float64 and Float32
for N in [Float64, Float32]
    # normalization
    H1 = Hyperplane(N[1e5, 2e5], N(3e5))
    H2 = normalize(H1)
    @test norm(H2.a) ≈ N(1) && H2.b == H1.b / norm(H1.a)
    @test normalize(H1, N(1)) == Hyperplane(N[1//3, 2//3], N(1))
    @test normalize(H1, N(Inf)) == Hyperplane(N[1//2, 1], N(3//2))

    # distance
    H = Hyperplane(N[1, -1], N(0))  # x <= y
    y = N[1, 1]  # closest point in the half-space
    for x in [N[2, 0], N[1, 1], N[0, 2]]
        @test distance(x, H) == distance(H, x) ≈ distance(x, y, p=N(2))
    end

    # sampling
    H = Hyperplane(N[1, -1], N(0))  # x = y
    for x in sample(H, 10)
        # membership
        @test x ∈ H
    end
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

    # tests that require Symbolics
    @static if isdefined(@__MODULE__, :Symbolics)
        vars = @variables x y
        @test Hyperplane(2x + 3y == 5) == Hyperplane([2.0, 3.0], 5.0)
        @test Hyperplane(2x + 3y == 5, N=Int) == Hyperplane([2, 3], 5)
        @test Hyperplane(2x + 3y == 5, vars) == Hyperplane([2.0, 3.0,], 5.0)
        @test Hyperplane(2x == 5y) == Hyperplane([2.0, -5.0,], 0.0)
        @test Hyperplane(2x == 5y, vars) == Hyperplane([2.0, -5.0,], 0.0)

        # test with sparse variables
        @variables x[1:5]
        @test Hyperplane(2x[1] + 5x[4] == 10., x) == Hyperplane([2.0, 0.0, 0.0, 5.0, 0.0], 10.0)

        # test passing a combination of operations
        vars = @variables x[1:2] t
        @test Hyperplane(x[1] == t, vars) == Hyperplane([1.0, 0.0, -1.0], 0.0)
    end
end
