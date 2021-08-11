for N in [Float64, Float32, Rational{Int}]
    c1 = N[1, 2]
    c2 = N[3, 4, 5]
    Z1 = Zonotope(c1, N[1 0; 0 1])
    Z2 = Zonotope(c2, N[1 0 0; 0 1 0; 0 0 1])
    H1 = Hyperrectangle(c1, N[1, 1])
    H2 = Hyperrectangle(c2, N[1, 1, 1])
    B1 = BallInf(c1, N(1))
    B2 = BallInf(c2, N(1))

    Z1Z2 = cartesian_product(Z1, Z2)
    Z1H2 = cartesian_product(Z1, H2)
    H1Z2 = cartesian_product(H1, Z2)
    H1H2 = cartesian_product(H1, H2)
    H1B2 = cartesian_product(H1, B2)
    B1B2 = cartesian_product(B1, B2)
    for X in [Z1Z2, Z1H2, H1Z2]
        @test X isa Zonotope{N}
    end
    for X in [H1H2, H1B2, B1B2]
        @test X isa Hyperrectangle{N}
    end
    if test_suite_polyhedra || N <: AbstractFloat
        for X in [Z1Z2, Z1H2, H1Z2, H1H2, H1B2]
            @test isequivalent(X, B1B2)
        end
    end

    # Universe is handled correctly (#2735)
    P = HalfSpace(N[1, 2], N(3))
    U = Universe{N}(2)
    @test isequivalent(cartesian_product(P, U), HalfSpace(N[1, 2, 0, 0], N(3)))
    @test isequivalent(cartesian_product(U, P), HalfSpace(N[0, 0, 1, 2], N(3)))
    @test cartesian_product(U, U) == Universe{N}(4)

    # hyperrectangle and Universe
    H = Hyperrectangle(N[1, 2], N[3, 4])
    HU = cartesian_product(H, U)
    @test HU isa HPolyhedron{N, SingleEntryVector{N}} && isequivalent(HU,
        HPolyhedron([HalfSpace(N[1, 0, 0, 0], N(4)),
                     HalfSpace(N[-1, 0, 0, 0], N(2)),
                     HalfSpace(N[0, 1, 0, 0], N(6)),
                     HalfSpace(N[0, -1, 0, 0], N(2))]))
    UH = cartesian_product(U, H)
    @test UH isa HPolyhedron{N, SingleEntryVector{N}} && isequivalent(UH,
        HPolyhedron([HalfSpace(N[0, 0, 1, 0], N(4)),
                     HalfSpace(N[0, 0, -1, 0], N(2)),
                     HalfSpace(N[0, 0, 0, 1], N(6)),
                     HalfSpace(N[0, 0, 0, -1], N(2))]))

    # HalfSpace and Universe
    H1 = HalfSpace(N[0, 4, 0], N(5))
    H2 = HalfSpace(sparsevec([2], N[4], 3), N(5))
    H3 = HalfSpace(SingleEntryVector(2, 3, N(4)), N(5))
    @test isequivalent(H1, H2) && isequivalent(H1, H3)
    H1a = cartesian_product(H1, Universe{N}(2))
    H1p = cartesian_product(Universe{N}(2), H1)
    @test isequivalent(H1a, cartesian_product(H2, Universe{N}(2))) &&
          isequivalent(H1a, cartesian_product(H3, Universe{N}(2))) &&
          isequivalent(H1a, HalfSpace(N[0, 4, 0, 0, 0], N(5)))
    @test isequivalent(H1p, cartesian_product(Universe{N}(2), H2)) &&
          isequivalent(H1p, cartesian_product(Universe{N}(2), H3)) &&
          isequivalent(H1p, HalfSpace(N[0, 0, 0, 4, 0], N(5)))

    # polyhedron and Universe
    P = HPolyhedron([HalfSpace(N[1], N(2))])
    PU = cartesian_product(P, U)
    @test isequivalent(PU, HPolyhedron([HalfSpace(N[1, 0, 0], N(2))]))
    UP = cartesian_product(U, P)
    @test isequivalent(UP, HPolyhedron([HalfSpace(N[0, 0, 1], N(2))]))

    # polytopes and polyhedra
    X = Interval(N(1), N(3))
    XX = BallInf(N[2, 2], N(1))
    P = convert(HPolytope, X)
    Q = convert(HPolyhedron, X)
    PP = cartesian_product(P, P)
    @test PP isa HPolytope{N} && isequivalent(PP, XX)
    QQ = cartesian_product(Q, Q)
    @test QQ isa HPolyhedron{N} && isequivalent(QQ, XX)
end
