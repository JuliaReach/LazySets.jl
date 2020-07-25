using LazySets: translate, ispermutation

for N in [Float64, Rational{Int}, Float32]
    # random line
    rand(Line)

    # construction
    l1 = Line(from=N[0, 1], to=N[1, 1]) # two points in the line
    l2 = Line(N[0, 1], N[1, 0]) # point and direction

    # the lines are the same modulo the sign of the normal vector
    @test l1.p ≈ l2.p && l1.d ≈ -l2.d

    # dimension
    @test dim(l1) == 2

    # support function
    @test ρ(N[0, 1], l1) == N(1)
    @test ρ(N[1, 0], l1) == N(Inf)

    # support vector
    @test σ(N[0, 1], l1) == N[0, 1]
    @test_throws ArgumentError σ(N[1, 0], l1) == N(Inf)

    # boundedness
    @test !isbounded(l1)

    # universality
    @test !isuniversal(l1)

    # isempty
    @test !isempty(l1)

    # an_element and membership
    an_element(l1) ∈ l1

    # constraints_list
    if N <: AbstractFloat
        @test ispermutation(constraints_list(l2),
                            [HalfSpace(N[0, 1], N(1)),    # y <= 1
                             HalfSpace(N[0, -1], N(-1))]) # y >= 1
    end

    # translation
    @test translate(l2, N[0, 1]) == Line(N[0, 2], N[1, 0])

    # distance
    distance(N[1, 0], l1) == N(1)
    distance(l1, N[1, 0]) == N(1)
    distance(Singleton(N[1, 0]), l1) == N(1)
    distance(l1, Singleton(N[1, 0])) == N(1)

    # concrete linear map
    if N <: AbstractFloat
        mirror = N[-1 0; 0 1]
        l = Line(from=N[0, 1], to=N[1, 1])
        @test isequivalent(linear_map(mirror, l), l)
        rot = N[0 -1; 1 0] # π/2 ccw rotation
        @test isequivalent(linear_map(rot, l), Line(N[-1, 0], N[0, 1]))
    end

    # constrained dimensions
    @test constrained_dimensions(l1) == [2]
end
