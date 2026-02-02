using LazySets, Test
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

for N in @tN([Float64, Float32, Rational{Int}])
    # interval \ half-space
    x = Interval(N(1), N(3))
    U = HalfSpace(N[1], N(2))
    @test difference(x, U) == Interval(N(2), N(3))
    U = HalfSpace(N[-1], N(-2))
    @test difference(x, U) == Interval(N(1), N(2))
    for U in (HalfSpace(N[1], N(3)), HalfSpace(N[-1], N(-1)))
        @test difference(x, U) == EmptySet{N}(1)
    end
    for U in (HalfSpace(N[1], N(1)), HalfSpace(N[-1], N(-3)))
        @test difference(x, U) == x
    end

    # polyhedron \ polyhedron
    P = VPolygon([N[-1, 0], N[1, 0], N[0, 2]])
    Q = VPolygon([N[0, 1], N[2, 2], N[1, 4]])
    R = difference(P, Q)
    # hard to test -> only sufficient test for result of implementation
    @test eltype(R) == N
    @test R == UnionSetArray([HPolytope(HalfSpace[HalfSpace(N[2, 1], N(2)),
                                                  HalfSpace(N[-2, 1], N(2)),
                                                  HalfSpace(N[0, -2], N(0)),
                                                  HalfSpace(N[3, -1], N(-1))]),
                              HPolytope(HalfSpace[HalfSpace(N[2, 1], N(2)),
                                                  HalfSpace(N[-1, 2], N(2)),
                                                  HalfSpace(N[-3, 1], N(1)),
                                                  HalfSpace(N[0, -2], N(0))])])
end
