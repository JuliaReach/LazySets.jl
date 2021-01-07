for N in [Float64, Rational{Int}, Float32]
    B = BallInf(ones(N, 2), N(3))
    H = Hyperrectangle(ones(N, 2), ones(N, 2))
    E = EmptySet{N}(2)

    # intersection of two sets
    I = Intersection(B, H)

    # swap
    I2 = swap(I)
    @test I.X == I2.Y && I.Y == I2.X

    # dim
    @test dim(I) == 2

    # support vector (currently throws an error)
    @test_throws ErrorException σ(ones(N, 2), I)

    # membership
    @test ones(N, 2) ∈ I && N[5, 5] ∉ I

    # boundedness (more tests in Float64-only section)
    @test isbounded(I)
    @test isbounded(Singleton(N[1]) ∩ HalfSpace(N[1], N(1)))

    # emptiness of intersection
    @test !isempty_known(I)
    @test !isempty(I)
    @test isempty_known(I)
    @test !isempty(I)

    # concretize
    @test concretize(I) == intersection(B, H)

    # conversion of intersection of polyhedral types to polyhedral types
    H1 = HalfSpace(N[1, 1], N(1))
    H2 = HalfSpace(N[-1, -1], N(1))
    P12 = convert(HPolyhedron, H1 ∩ H2)
    @test P12 == HPolyhedron([H1, H2])

    # simplification of lazy intersection of half-spaces to polyhedron
    @test H1 ∩ H2 == HPolyhedron([H1, H2])
    H3 = HalfSpace(N[1, 0], N(3))
    @test H1 ∩ H2 ∩ H3 == HPolyhedron([H1, H2, H3])
    # make lazy intersection of half-spaces by specifying the cache
    Ih12 = Intersection(H1, H2, cache=LazySets.IntersectionCache())
    @test Ih12.X == H1 && Ih12.Y == H2

    # vertices
    I2 = Intersection(BallInf(zeros(N, 2), N(2)), BallInf(3 * ones(N, 2), N(2)))
    @test ispermutation(vertices_list(I2),
        vertices_list(BallInf(fill(N(3//2) , 2), N(1//2))))

    # =================
    # IntersectionArray
    # =================

    # relation to base type (internal helper functions)
    @test LazySets.array_constructor(Intersection) == IntersectionArray
    @test LazySets.is_array_constructor(IntersectionArray)

    # intersection of an array of sets
    IArr = IntersectionArray([B, H])

    # dim
    @test dim(IArr) == 2

    # support vector (currently throws an error)
    @test_throws ErrorException σ(ones(N, 2), IArr)

    # boundedness
    @test isbounded(IArr)
    @test isbounded(IntersectionArray([Singleton(N[1]), HalfSpace(N[1], N(1))]))
    # the following tests crash because ρ(::IntersectionArray) is not implemented yet
    @test_throws ErrorException isbounded(IntersectionArray([HalfSpace(N[1], N(1)), HalfSpace(N[1], N(-1))]))
    @test_throws ErrorException !isbounded(IntersectionArray([HalfSpace(ones(N, 2), N(1)), HalfSpace(ones(N, 2), N(-1))]))

    # isempty
    @test_throws MethodError isempty(IArr)

    # membership
    @test ones(N, 2) ∈ IArr && N[5, 5] ∉ IArr

    # array getter
    v = Vector{LazySet{N}}()
    @test array(IntersectionArray(v)) ≡ v

    # constructor with size hint and type
    IntersectionArray(10, N)

    # concretize
    @test concretize(IArr) == intersection(B, H)

    # conversion of intersection of polyhedral types to polyhedral types
    P12 = convert(HPolyhedron, IntersectionArray([H1, H2]))
    @test P12 == HPolyhedron([H1, H2])

    # ================
    # common functions
    # ================

    # absorbing element
    @test absorbing(Intersection) == absorbing(IntersectionArray) == EmptySet
    @test I ∩ E == E ∩ I == IArr ∩ E == E ∩ IArr == E ∩ E == E
end

# ======================
# Tests for Float64 only
# ======================
for N in [Float64]
    # constraints_list for polytopic intersection
    B = BallInf(ones(N, 2), N(3))
    H = Hyperrectangle(ones(N, 2), ones(N, 2))
    I = B ∩ H
    IArr = IntersectionArray([B, H])
    clist1 = constraints_list(I)
    clist2 = constraints_list(IArr)
    @test ispermutation(clist1, clist2) &&
          ispermutation(clist1, [HalfSpace(N[1, 0], N(2)),
                                 HalfSpace(N[0, 1], N(2)),
                                 HalfSpace(N[-1, 0], N(0)),
                                 HalfSpace(N[0, -1], N(0))])

    # HalfSpace vs. Ball1 intersection
    X = Ball1(zeros(2), N(1));
    d = normalize(N[1, 0])

    # flat intersection at x = 1
    H = HalfSpace(N[-1, 0], N(-1)); # x >= 1

    # default algorithm
    @test ρ(d, X ∩ H) == ρ(d, H ∩ X) == N(1)

    # intersection at x = 0
    H = HalfSpace(N[1, 0], N(0)); # x <= 0

    # default algorithm
    @test ρ(d, X ∩ H) < 1e-6 && ρ(d, H ∩ X) < 1e-6

    # specify  line search algorithm
    @test ρ(d, X ∩ H, algorithm="line_search") < 1e-6 &&
        ρ(d, H ∩ X, algorithm="line_search") < 1e-6

    # boundedness
    @test isbounded(HalfSpace(N[1], N(1)) ∩ HalfSpace(N[-1], N(1)))
    @test !isbounded(HalfSpace(ones(N, 2), N(1)) ∩ HalfSpace(-ones(N, 2), N(1)))

    # HalfSpace vs. Ball2 intersection
    B2 = Ball2(zeros(2), N(1));
    @test ρ(d, B2 ∩ H) < 1e-6 && ρ(d, H ∩ B2) < 1e-6

    # Ball1 vs. Hyperplane intersection
    H = Hyperplane(N[1, 0], N(0.5)); # x = 0.5
    @test isapprox(ρ(d, X ∩ H, algorithm="line_search"), N(0.5), atol=1e-6)
    # For the projection algorithm, if the linear map is taken lazily we can use Ball1
    @test isapprox(ρ(d, X ∩ H, algorithm="projection", lazy_linear_map=true), N(0.5), atol=1e-6)
    # But the default is to take the linear map concretely; in this case, we *may*
    # need Polyhedra (in the general case), for the concrete linear map. As a valid workaround
    # if we don't want to load Polyhedra here is to convert the given set to a polygon in V-representation
    @test isapprox(ρ(d, convert(VPolygon, X) ∩ H, algorithm="projection", lazy_linear_map=false), N(0.5), atol=1e-6)

    # =====================
    # concrete operations
    # =====================
    cap =  HPolytope([HalfSpace(N[1], N(1))]) ∩ HPolytope([HalfSpace(N[-1], N(1))])  # x <= 1 && x >= -1
    p = linear_map(reshape([N(1/2)], 1, 1), cap)
    @test N[-0.5] ∈ p && N[0.5] ∈ p && (N[1.0] ∉ p || N[1.0] ∉ p)
end
