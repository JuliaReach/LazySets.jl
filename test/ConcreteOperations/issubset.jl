using LazySets, Test
IA = LazySets.IA
using LazySets.IA: IntervalBox
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

for N in @tN([Float64, Float32, Rational{Int}])
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

    # AbstractZonotope in AbstractHyperrectangle
    Z = Zonotope(N[0, 0], N[1 1; -1 1])
    H = Hyperrectangle(; low=N[-2, -2], high=N[2, 2])
    @test Z ⊆ H
    res, w = ⊆(Z, H, true)
    @test res && w isa Vector{N} && isempty(w)
    H = Hyperrectangle(; low=N[-2, -2], high=N[2, 0])
    @test !(Z ⊆ H)
    @test_broken ⊆(Z, H, true) isa Tuple  # TODO this should work

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
    @test_throws ArgumentError issubset(Pnc, P, true)  # not implemented; parser requires `issubset`

    # zonotope in polyhedron
    Z = Zonotope(N[1, 1], N[1 0; 0 1])
    P = HPolyhedron([HalfSpace(N[1, 0], N(3)), HalfSpace(N[-1, 0], N(0)),
                     HalfSpace(N[0, 1], N(3)), HalfSpace(N[0, -1], N(0))])
    res, w = ⊆(Z, P, true)
    @test Z ⊆ P && res && w == N[]
    P = HPolyhedron([HalfSpace(N[1, 0], N(1)), HalfSpace(N[-1, 0], N(0)),
                     HalfSpace(N[0, 1], N(3)), HalfSpace(N[0, -1], N(0))])
    @test !(Z ⊆ P)
    @test_throws ArgumentError issubset(Z, P, true)  # not implemented; parser requires `issubset`
    # corner case: no generator
    Z = Zonotope(N[2], zeros(N, 1, 0))
    P = convert(HPolytope, Interval(N(1), N(3)))
    @test Z ⊆ P
end

for N in [Float64]
    # rounding error
    Z = Zonotope([-8.0, 8.0],
                 [-1.0000000000000002 0.7000000000000002;
                  1.0000000000000002 -0.7000000000000002])
    L = LineSegment([-9.7, 9.7], [-6.3, 6.3])
    @test Z ⊆ L
end
