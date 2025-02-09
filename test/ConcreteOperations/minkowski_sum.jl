for N in [Float64, Float32, Rational{Int}]
    X = Interval(N(1), N(2))
    Y = X + X
    if test_suite_polyhedra
        minkowski_sum(X, Y)
        minkowski_sum(X, Y; algorithm=Polyhedra.FourierMotzkin())
    end

    for B in (Ball1, BallInf)
        B1 = B(N[1, 2], N(3))
        B2 = B(N[4, 5], N(6))
        @test minkowski_sum(B1, B2) == B(N[5, 7], N(9))
    end

    H1 = Hyperrectangle(N[1, 2], N[3, 4])
    H2 = Hyperrectangle(N[5, 6], N[7, 8])
    @test minkowski_sum(H1, H2) == Hyperrectangle(N[6, 8], N[10, 12])
    PZ = minkowski_sum(convert(SparsePolynomialZonotope, H1), H2)
    # equality is not required but approximates the equivalence check
    @test PZ == SparsePolynomialZonotope(N[6, 8], N[3 0; 0 4], N[7 0; 0 8], [1 0; 0 1], 1:2)
    PZ = minkowski_sum(H1, convert(SparsePolynomialZonotope, H2))
    @test PZ == SparsePolynomialZonotope(N[6, 8], N[7 0; 0 8], N[3 0; 0 4], [1 0; 0 1], 1:2)

    # EmptySet / Universe
    E = EmptySet{N}(2)
    U = Universe{N}(2)
    Z = ZeroSet{N}(2)
    B = Ball1(N[0, 0], N(1))
    for E2 in (minkowski_sum(E, U), minkowski_sum(U, E), minkowski_sum(E, Z), minkowski_sum(Z, E),
               minkowski_sum(E, B), minkowski_sum(B, E))
        @test E2 isa EmptySet{N} && E2 == E
    end
    for U2 in (minkowski_sum(U, Z), minkowski_sum(Z, U), minkowski_sum(U, B), minkowski_sum(B, U))
        @test U2 isa Universe{N} && U2 == U
    end
end

for N in [Float64, Float32]
    B1 = Ball2(N[1, 2], N(3))
    B2 = Ball2(N[4, 5], N(6))
    @test minkowski_sum(B1, B2) == Ball2(N[5, 7], N(9))

    B1 = Ballp(N(3), N[1, 2], N(3))
    B2 = Ballp(N(3), N[4, 5], N(6))
    @test minkowski_sum(B1, B2) == Ballp(N(3), N[5, 7], N(9))
end
