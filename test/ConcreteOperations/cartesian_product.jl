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

    # Universe with hyperrectangle (#2741)
    H = Hyperrectangle(N[1, 2], N[3, 4])
    HU = cartesian_product(H, U)
    @test HU isa HPolyhedron{N, SingleEntryVector{N}} && isequivalent(HU,
        HPolyhedron([HalfSpace(N[1, 0, 0, 0], N(4)),
                     HalfSpace(N[-1, 0, 0, 0], N(2)),
                     HalfSpace(N[0, 1, 0, 0], N(6)),
                     HalfSpace(N[0, -1, 0, 0], N(2))]))
end
