using Optim

for N in [Float64, Rational{Int}, Float32]
    B = BallInf(ones(N, 2), N(3))
    H = Hyperrectangle(ones(N, 2), ones(N, 2))
    E = EmptySet{N}()

    # intersection of two sets
    I = Intersection(B, H)

    # dim
    @test dim(I) == 2

    # support vector (currently throws an error)
    @test_throws ErrorException σ(ones(N, 2), I)

    # membership
    @test ∈(ones(N, 2), I) && !∈(N[5, 5], I)

    # boundedness (more tests in Float64-only section)
    @test isbounded(I)
    @test isbounded(Singleton(N[1]) ∩ HalfSpace(N[1], N(1)))

    # emptiness of intersection
    @test !isempty_known(I)
    @test !isempty(I)
    @test isempty_known(I)
    @test !isempty(I)

    # constraints_list for polytopic intersection
    @test ispermutation(constraints_list(I),
                        [HalfSpace{Float64}(N[1, 0], N(2)),
                         HalfSpace{Float64}(N[0, 1], N(2)),
                         HalfSpace{Float64}(N[-1, 0], N(0)),
                         HalfSpace{Float64}(N[0, -1], N(0))])

    # =================
    # IntersectionArray
    # =================

    # relation to base type (internal helper functions)
    @test LazySets.array_constructor(Intersection) == IntersectionArray
    @test LazySets.is_array_constructor(IntersectionArray)

    # intersection of an array of sets
    IA = IntersectionArray([B, H])

    # dim
    @test dim(IA) == 2

    # support vector (currently throws an error)
    @test_throws ErrorException σ(ones(N, 2), IA)

    # boundedness
    @test isbounded(IA)
    @test isbounded(IntersectionArray([Singleton(N[1]), HalfSpace(N[1], N(1))]))
    # the following tests crash because ρ(::IntersectionArray) is not implemented yet
    @test_throws ErrorException isbounded(IntersectionArray([HalfSpace(N[1], N(1)), HalfSpace(N[1], N(-1))]))
    @test_throws ErrorException !isbounded(IntersectionArray([HalfSpace(ones(N, 2), N(1)), HalfSpace(ones(N, 2), N(-1))]))

    # isempty
    @test_throws MethodError isempty(IA)

    # membership
    @test ∈(ones(N, 2), IA) && !∈(N[5, 5], IA)

    # array getter
    v = Vector{LazySet{N}}()
    @test array(IntersectionArray(v)) ≡ v

    # constructor with size hint and type
    IntersectionArray(10, N)

    # constraints_list for polytopic intersection
    @test ispermutation(constraints_list(IA),
                        [HalfSpace{Float64}(N[1, 0], N(2)),
                         HalfSpace{Float64}(N[0, 1], N(2)),
                         HalfSpace{Float64}(N[-1, 0], N(0)),
                         HalfSpace{Float64}(N[0, -1], N(0))])

    # ================
    # common functions
    # ================

    # absorbing element
    @test absorbing(Intersection) == absorbing(IntersectionArray) == EmptySet
    @test I ∩ E == E ∩ I == IA ∩ E == E ∩ IA == E ∩ E == E
end

# ======================
# Tests for Float64 only
# ======================
for N in [Float64]
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
end
