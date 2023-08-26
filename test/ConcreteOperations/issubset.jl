for N in [Float64, Float32, Rational{Int}]
    # interval in union of intervals
    X = Interval(N(0), N(1))

    # included
    YU = Interval(N(1 // 2), N(2)) ∪ Interval(N(0), N(2 // 3))
    YU_hr = convert(Hyperrectangle, YU.X) ∪ convert(Hyperrectangle, YU.Y)
    Yarr = UnionSetArray([YU.X, YU.Y])
    for Y in [YU, YU_hr, Yarr]
        res, w = ⊆(X, Y, true)
        @test X ⊆ Y && res && w == N[]
    end

    # not included
    YU = Interval(N(1 // 2), N(5 // 6)) ∪ Interval(N(0), N(1 // 3))
    YU_hr = convert(Hyperrectangle, YU.X) ∪ convert(Hyperrectangle, YU.Y)
    Yarr = UnionSetArray([YU.X, YU.Y])
    for Y in [YU, YU_hr, Yarr]
        res, w = ⊆(X, Y, true)
        @test !(X ⊆ Y) && !res && w ∈ X && w ∉ Y
    end

    # using IA types
    X = interval(N(0), N(1)) # IA
    Y = Interval(N(-1), N(2))
    @test X ⊆ Y
    @test !(Y ⊆ X)

    X = IntervalBox(interval(N(0), N(1)), interval(N(0), N(1)))
    Y = Hyperrectangle(; low=[N(-1), N(-1)], high=[N(2), N(2)])
    @test X ⊆ Y
    @test !(Y ⊆ X)

    # interval with union of unbounded sets
    X = Interval(N(1), N(3))
    for H in (HalfSpace(N[1], N(2)), HalfSpace(N[-1], N(-2)))
        U = UnionSetArray([H])
        res, w = ⊆(X, U, true)
        @test !(X ⊆ U) && !res && w ∈ X && w ∉ U
    end
    X = Interval(N(1), N(2))
    for H in (HalfSpace(N[1], N(2)), HalfSpace(N[-1], N(-1)))
        U = UnionSetArray([H])
        res, w = ⊆(X, U, true)
        @test X ⊆ U && res && w == N[]
    end

    B = Hyperrectangle(N[0, 0], N[1, 1])
    U = UnionSetArray([Hyperrectangle(N[0, 1], N[3, 1]), Hyperrectangle(N[0, -1], N[3, 1])])
    res, w = ⊆(B, U, true)
    @test B ⊆ U && res && w == N[]
    U = UnionSetArray([Hyperrectangle(N[0, 1], N[3, 1]), Hyperrectangle(N[3, -1], N[3, 1])])
    res, w = ⊆(B, U, true)
    @test !(B ⊆ U) && !res && w ∈ B && w ∉ U

    # nonconvex set in polyhedron
    Pnc = Polygon([N[0, 0], N[0, 2], N[2, 2], N[2, 0], N[1, 1]])
    P = HPolyhedron([HalfSpace(N[1, 0], N(3)), HalfSpace(N[-1, 0], N(0)),
                     HalfSpace(N[0, 1], N(3)), HalfSpace(N[0, -1], N(0))])
    res, w = ⊆(Pnc, P, true)
    @test Pnc ⊆ P && res && w == N[]
    P = HPolyhedron([HalfSpace(N[1, 0], N(1)), HalfSpace(N[-1, 0], N(0)),
                     HalfSpace(N[0, 1], N(1)), HalfSpace(N[0, -1], N(0))])
    @test !(Pnc ⊆ P)
    @test_throws ArgumentError issubset(Pnc, P, true)  # not implemented

    # zonotope in polyhedron
    Z = Zonotope(N[1, 1], N[1 0; 0 1])
    P = HPolyhedron([HalfSpace(N[1, 0], N(3)), HalfSpace(N[-1, 0], N(0)),
                     HalfSpace(N[0, 1], N(3)), HalfSpace(N[0, -1], N(0))])
    res, w = ⊆(Z, P, true)
    @test Z ⊆ P && res && w == N[]
    P = HPolyhedron([HalfSpace(N[1, 0], N(1)), HalfSpace(N[-1, 0], N(0)),
                     HalfSpace(N[0, 1], N(3)), HalfSpace(N[0, -1], N(0))])
    @test !(Z ⊆ P)
    @test_throws ArgumentError issubset(Z, P, true)  # not implemented
end
