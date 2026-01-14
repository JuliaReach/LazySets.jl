using LazySets, Test
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

for N in @tN([Float64, Float32, Rational{Int}])
    # Interval and HalfSpace
    X = Interval(N(1), N(2))
    H = HalfSpace(N[1], N(1.5))
    @test intersection(X, H) == Interval(N(1), N(1.5))
    H = HalfSpace(N[-2], N(-5))
    @test intersection(X, H) == EmptySet{N}(1)
    H = HalfSpace(N[2], N(5))
    @test intersection(X, H) == X
    # Interval and Hyperplane
    H = Hyperplane(N[2], N(3))
    @test intersection(X, H) == Singleton(N[1.5])
    H = Hyperplane(N[-1], N(-3))
    @test intersection(X, H) == EmptySet{N}(1)
    # Interval and Ball1
    Y = Ball1(N[2], N(0.5))
    @test intersection(X, Y) == Interval(N(1.5), N(2))
    # Interval and ConvexHull
    Y = ConvexHull(Singleton(N[-5]), Singleton(N[-1]))
    @test intersection(X, Y) == EmptySet{N}(1)
    # Interval and UnionSetArray
    U = UnionSetArray([Interval(N[-1, 2])])
    @test intersection(U, X) == intersection(X, U) == X

    # AbstractSingleton and AbstractHyperrectangle
    S = Singleton(N[2, -1])
    H = Hyperrectangle(N[2, 0], N[1, 1])
    X = intersection(S, H)
    @test X isa LazySet{N} && isequivalent(X, S)
    H = Hyperrectangle(N[4, 0], N[1, 1])
    X = intersection(S, H)
    @test X isa LazySet{N} && isequivalent(X, EmptySet{N}(2))
    # AbstractSingleton and LazySet
    P = Polygon([N[2, -1]])
    X = intersection(S, P)
    @test X isa LazySet{N} && isequivalent(X, S)
    P = Polygon([N[3, 1]])
    X = intersection(S, P)
    @test X isa LazySet{N} && isequivalent(X, EmptySet{N}(2))

    # `_intersection_poly`
    # test specializations in 1D and 2D plus the general case in 3D
    # TODO add tests for keyword arguments
    for n in 1:3
        a = zeros(N, n)
        a[1] = 1
        H1 = HalfSpace(a, N(1))  # x <= 1
        H4 = HalfSpace(a, N(0))  # x <= 0
        a = zeros(N, n)
        a[1] = -1
        H2 = HalfSpace(a, N(0))  # x >= 0
        H3 = HalfSpace(a, N(-2))  # x >= 2
        B = convert(HPolytope, BallInf(zeros(N, n), N(1)))

        # - unbounded inputs
        #   - bounded output
        X = intersection(H1, H2)
        @test isequivalent(X, HPolyhedron([H1, H2]))
        #   - empty output
        X = intersection(H1, H3)
        @test X == EmptySet{N}(n)
        #   - unbounded output
        X = intersection(H1, H4)
        @test isequivalent(X, H4)
        # - bounded inputs
        #   - bounded output
        X = intersection(H1, B)
        @test isequivalent(X, B)
        #   - empty output
        X = intersection(H3, B)
        @test X == EmptySet{N}(n)
    end
    # 1D bounded from below only
    H1 = HalfSpace(N[-1], N(0))  # x >= 0
    H2 = HalfSpace(N[-1], N(-1))  # x >= 1
    X = intersection(H1, H2)
    @test isequivalent(X, H2)
end

for N in @tN([Float64, Rational{Int}])
    # Zonotope and HalfSpace
    Z = Zonotope(N[1, 2], N[1 3; 2 4])
    H = HalfSpace(N[1, 0], N(-4))  # disjoint
    X = intersection(H, Z)
    @test X isa EmptySet{N} && X == EmptySet{N}(2)
    H = HalfSpace(N[1, 0], N(6))  # fully contained
    X = intersection(H, Z)
    @test X == Z
    H = HalfSpace(N[1, 0], N(0))  # overlapping
    X = VPolygon([N[-3, -4], N[0, 0], N[0, 4 // 3], N[-1, 0]])
    @test isequivalent(intersection(H, Z), X)
end
