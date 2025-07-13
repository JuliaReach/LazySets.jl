for N in @tN([Float64, Float32, Rational{Int}])
    # intervals, nonempty difference
    X = Interval(N(1), N(3))
    Y = Interval(N(2), N(3))
    @test minkowski_difference(X, Y) == Interval(N(-1), N(0))

    # intervals, empty difference
    X = Interval(N(1), N(3))
    Y = Interval(N(2), N(5))
    D = minkowski_difference(X, Y)
    @test D isa EmptySet{N} && dim(D) == 1

    # hyperrectangles, nonempty difference
    H1 = Hyperrectangle(N[2, 3], N[1, 3])
    H2 = Hyperrectangle(N[5 / 2, -1], N[1 / 2, 2])
    @test minkowski_difference(H1, H2) == Hyperrectangle(N[-1 / 2, 4], N[1 / 2, 1])

    # hyperrectangles, empty difference
    H1 = Hyperrectangle(N[2, 3], N[1, 3])
    H2 = Hyperrectangle(N[5, -1], N[2, 2])
    D = minkowski_difference(H1, H2)
    @test D isa EmptySet{N} && dim(D) == 2

    # zonotopes in 2D (taken from [Althoff15; Fig. 2 in the 2022 version](@citet))
    Zm = Zonotope(N[1, 1], N[1 0 1; 0 1 1])

    if N isa AbstractFloat
        Zs1 = Zonotope(N[0, 0], N[1/2 0; -1/5 1/5])
        D1 = minkowski_difference(Zm, Zs1)
        @test D1 ≈ Zonotope(N[1, 1], N[1/2 0 1; 0 3/5 1])
    end

    Zs2 = Zonotope(N[0, 0], N[1/2 0; -1/2 1/2])
    D2 = minkowski_difference(Zm, Zs2)
    @test D2 ≈ Zonotope(N[1, 1], N[1/2 1; 0 1])

    Zs3 = Zonotope(N[0, 0], N[2 0; -1/2 1/2])
    D3 = minkowski_difference(Zm, Zs3)
    @test D3 isa EmptySet{N} && dim(D3) == 2

    # zonotope in 3D
    Zm = Zonotope(N[1, 1, 1], N[1 0 0; 0 1 0; 0 0 1])
    Zs = Zonotope(N[0, 0, 1], hcat(N[1 / 2; 0; 0]))
    D = minkowski_difference(Zm, Zs)
    @test isequivalent(D,
                       HPolytope([HalfSpace(N[0, 0, 1], N(1)), HalfSpace(N[0, 0, -1], N(1)),
                                  HalfSpace(N[0, 1, 0], N(2)), HalfSpace(N[0, -1, 0], N(0)),
                                  HalfSpace(N[1, 0, 0], N(3 / 2)),
                                  HalfSpace(N[-1, 0, 0], N(-1 / 2))]))
    if N == Float64
        @test isequivalent(D, Zonotope(N[1, 1, 0], N[1/2 0 0; 0 1 0; 0 0 1]))
    end
end

for N in [Float64]
    if test_suite_polyhedra
        # concrete Minkowski difference for unbounded P (HPolyhedron)
        mx1 = N(2)
        mx2 = N(5)
        P3 = HPolyhedron(N[mx1 0; 0 mx2], N[3, 3])
        radius = 2
        Q3 = Ball1(N[0, 0], N(radius))
        C3 = minkowski_difference(P3, Q3)
        C3_res = HPolyhedron(N[mx1 0; 0 mx2], N[3 - mx1 * radius, 3 - mx2 * radius])
        @test C3 ⊆ C3_res && C3_res ⊆ C3
    end
end
