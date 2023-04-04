for N in [Float64, Rational{Int}, Float32]
    # random empty set
    rand(EmptySet)

    E = EmptySet{N}(2)
    B = BallInf(ones(N, 2), N(1))

    # dim
    @test dim(E) == 2

    # support vector
    @test_throws ErrorException σ(N[0], E)

    # boundedness
    @test isbounded(E) && isboundedtype(typeof(E))

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

    # concrete rectification
    @test rectify(E) == E
end

# tests that only work with Float64 and Float32
for N in [Float64, Float32]
    B = Ball2(N[0, 0], N(1))
    E = EmptySet{N}(2)

    @test isdisjoint(E, B) && isdisjoint(B, E)
end

# default Float64 constructor
@test EmptySet(2) == ∅(2) == EmptySet{Float64}(2)

# intersection
for X in filter(!isoperationtype, LazySets.subtypes(LazySet, true))
    if X ∈ [RotatedHyperrectangle, Star, Polygon, DensePolynomialZonotope,
            SparsePolynomialZonotope, SimpleSparsePolynomialZonotope]
        # missing rand()
        continue
    end
    Y = rand(X)
    E = EmptySet(dim(Y))
    @test intersection(Y, E) == intersection(E, Y) == E
end
