for N in [Float64, Rational{Int}, Float32]
    # random line
    rand(Line2D)

    # construction
    a1 = N[0, 1]
    b1 = N(1)
    l1 = Line2D(a1, b1)
    l2 = Line2D(N[1, 0], N(1))
    l3 = Line2D(N[0, 1], N(2))
    l4 = Line2D(N[1, 1], N(0))
    l1_copy = Line2D(N[0, 1], N(1))

    # alternative construction from two points
    # line with positive slope
    L = Line2D([1.0, 0.0], [2.0, 0.5])
    @test [1.0, 0.0] ∈ L &&
          [2.0, 0.5] ∈ L &&
          [2.0, 1.5] ∉ L

    # vertical line
    L = Line2D([1.0, 0.0], [1.0, 2.5])
    @test [1.0, 0.0] ∈ L &&
          [1.0, 2.5] ∈ L &&
          [1.0, 10.0] ∈ L
          [10.0, 10.0] ∉ L

    # corner case: zero normal vector
    @test_throws AssertionError Line2D(N[0, 0], N(1))

    # dimension
    @test dim(l1) == 2

    # support vector
    σ(N[0, 1], l1)
    σ(N[1, 0], l2)
    σ(N[0, 1], l3)

    # boundedness
    @test !isbounded(l1)

    # universality
    @test !isuniversal(l1)
    res, w = isuniversal(l1, true)
    @test !res && w ∉ l1

    # isempty
    @test !isempty(l1)

    # an_element and membership
    an_element(l1) ∈ l1
    an_element(l2) ∈ l2
    an_element(l3) ∈ l3
    an_element(l4) ∈ l4

    # constrained dimensions
    @test constrained_dimensions(l1) == [2]
    @test constrained_dimensions(l4) == [1, 2]

    # constraints_list
    @test ispermutation(constraints_list(l1),
                        [HalfSpace(a1, b1), HalfSpace(-a1, -b1)])

    # concrete intersection
    cap11 = intersection(l1, l1_copy)
    cap12 = intersection(l1, l2)
    cap13 = intersection(l1, l3)
    @test cap11 isa Line2D && cap11.a == l1.a && cap11.b == l1.b
    @test cap12 isa Singleton && element(cap12) == N[1, 1]
    @test cap13 isa EmptySet{N}

    # concrete linear map of a line
    L = Line2D(N[1, -1], N(0)) # x = y
    M = N[1 0; 0 0] # non-invertible matrix
    # projection is y = 0
    if test_suite_polyhedra
        lm = linear_map(M, L)
        if N == Float32 || N == Float64
            @test lm isa Line2D{Float64}
            @test lm.a ≈ N[0, -1] && lm.b ≈ N(0)

            # returned set is universal
            @test linear_map(N[1 1], L) == Universe{N}(1)
        elseif N == Rational{Int}
            @test lm isa Line2D{Rational{BigInt}}
            @test lm.a == N[0//1, -1//1] && lm.b == N(0//1)
        end
    end
    M = N[2 2; 0 1] # invertible matrix
    lm = linear_map(M, L)
    @test lm isa Line2D{N, Vector{N}}
    @test lm.a ≈ N[1/2, -2] && lm.b ≈ N(0)

    # translation
    @test translate(l1, N[1, 2]) == Line2D(a1, N(3))
end
