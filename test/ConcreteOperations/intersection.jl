using LazySets, Test
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

for N in @tN([Float64, Float32, Rational{Int}])
    # Interval with HalfSpace
    X = Interval(N(1), N(2))
    H = HalfSpace(N[1], N(1.5))
    @test intersection(X, H) == Interval(N(1), N(1.5))
    H = HalfSpace(N[-2], N(-5))
    @test intersection(X, H) == EmptySet{N}(1)
    H = HalfSpace(N[2], N(5))
    @test intersection(X, H) == X
    # Interval with Hyperplane
    H = Hyperplane(N[2], N(3))
    @test intersection(X, H) == Singleton(N[1.5])
    H = Hyperplane(N[-1], N(-3))
    @test intersection(X, H) == EmptySet{N}(1)
    # Interval with Ball1
    Y = Ball1(N[2], N(0.5))
    @test intersection(X, Y) == Interval(N(1.5), N(2))
    # Interval with ConvexHull
    Y = ConvexHull(Singleton(N[-5]), Singleton(N[-1]))
    @test intersection(X, Y) == EmptySet{N}(1)
    # Interval with UnionSetArray
    U = UnionSetArray([Interval(N[-1, 2])])
    @test intersection(U, X) == intersection(X, U) == X
end

for N in @tN([Float64, Rational{Int}])
    # Zonotope with HalfSpace
    Z = Zonotope(N[1, 2], N[1 3; 2 4])
    H = HalfSpace(N[1, 0], N(-4))  # disjoint
    X = intersection(H, Z)
    @test X isa EmptySet{N} && X == EmptySet{N}(2)
    H = HalfSpace(N[1, 0], N(6))  # fully contained
    X = intersection(H, Z)
    @test isidentical(X, Z)
    H = HalfSpace(N[1, 0], N(0))  # overlapping
    X = VPolygon([N[-3, -4], N[0, 0], N[0, 4 // 3], N[-1, 0]])
    @test isequivalent(intersection(H, Z), X)
end
