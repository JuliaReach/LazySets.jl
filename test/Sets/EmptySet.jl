for N in [Float64, Rational{Int}, Float32]
    # random empty set
    rand(EmptySet)

    E = EmptySet{N}(2)
    B = BallInf(ones(N, 2), N(1))

    # dim
    @test dim(E) == 2

    # support function & support vector
    @test_throws AssertionError ρ(N[0], E)
    @test_throws AssertionError σ(N[0], E)

    # boundedness
    @test isbounded(E) && isboundedtype(typeof(E))

    # isconvextype
    @test isconvextype(typeof(E))

    # isuniversal
    res, w = isuniversal(E, true)
    @test !isuniversal(E) && !res && w ∉ E

    # membership
    @test_throws AssertionError N[0] ∈ E
    @test N[0, 0] ∉ E

    # subset
    @test E ⊆ B && ⊆(E, B, true)[1]
    subset, point = ⊆(B, E, true)
    @test B ⊈ E && !subset && point ∈ B
    @test E ⊆ E && ⊆(E, E, true)[1]

    # emptiness check
    @test isempty(E)

    # an_element/norm/radius/diameter functions
    @test_throws ArgumentError an_element(E)
    @test_throws ArgumentError norm(E)
    @test_throws ArgumentError radius(E)
    @test_throws ArgumentError diameter(E)

    # vertices / vertices_list
    @test collect(vertices(E)) == vertices_list(E) == Vector{Vector{N}}()
    @test typeof(collect(vertices(E))) == typeof(vertices_list(E)) ==
          Vector{Vector{N}}

    # linear map of an empty set
    @test linear_map(ones(N, 2, 2), E) == E
    @test linear_map(ones(N, 3, 2), E) == EmptySet{N}(3)

    # translation
    @test translate(E, N[1, 2]) == E
    @test translate!(E, N[1, 2]) == E

    # disjointness
    for X in [B, Singleton(N[0, 0])]
        @test isdisjoint(E, X) && isdisjoint(X, E)
    end

    # projection
    @test project(EmptySet{N}(5), [1, 4, 5]) == EmptySet{N}(3)

    # volume
    @test volume(E) == zero(N)

    # complement
    C = complement(E)
    @test C == Universe{N}(2) && C isa Universe{N}

    # concrete rectification
    @test rectify(E) == E

    # reflect
    @test reflect(E) == E

    # scale/scale!
    E2 = copy(E)
    scale!(N(2), E2)
    @test scale(N(2), E) == E2 == EmptySet{N}(dim(E))

    # area
    @test area(E) == N(0)

    # chebyshev_center_radius
    @test_throws ArgumentError chebyshev_center_radius(E)

    # low/high
    @test_throws ArgumentError low(E)
    @test_throws ArgumentError low(E, 1)
    @test_throws ArgumentError high(E)
    @test_throws ArgumentError high(E, 1)

    # isdisjoint
    @test isdisjoint(E, B) && isdisjoint(B, E) && isdisjoint(E, E)
    for (res, w) in (isdisjoint(E, B, true), isdisjoint(B, E, true), isdisjoint(E, E, true))
        @test res && w == N[]
    end

    # convex_hull
    @test convex_hull(E, E) == E
end

# default Float64 constructor
@test EmptySet(2) == ∅(2) == EmptySet{Float64}(2)

# intersection
for X in filter(!isoperationtype, LazySets.subtypes(LazySet, true))
    if X ∈ [Star, Polygon, DensePolynomialZonotope,
            SparsePolynomialZonotope, SimpleSparsePolynomialZonotope]
        # missing rand()
        continue
    end
    Y = rand(X)
    E = EmptySet(dim(Y))
    @test intersection(Y, E) == intersection(E, Y) == E
end
