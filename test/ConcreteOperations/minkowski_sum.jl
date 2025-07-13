for N in @tN([Float64, Float32, Rational{Int}])
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
    SPZ = convert(SparsePolynomialZonotope, H1)
    for P in (SPZ, convert(SSPZ, SPZ))
        P = minkowski_sum(P, H2)
        # equality is not required but approximates the equivalence check
        @test P == SparsePolynomialZonotope(N[6, 8], N[3 0; 0 4], N[7 0; 0 8], [1 0; 0 1], 1:2)
        P = minkowski_sum(H1, convert(SparsePolynomialZonotope, H2))
        @test P == SparsePolynomialZonotope(N[6, 8], N[7 0; 0 8], N[3 0; 0 4], [1 0; 0 1], 1:2)
    end
end

for N in @tN([Float64, Float32])
    B1 = Ball2(N[1, 2], N(3))
    B2 = Ball2(N[4, 5], N(6))
    @test minkowski_sum(B1, B2) == Ball2(N[5, 7], N(9))

    B1 = Ballp(N(3), N[1, 2], N(3))
    B2 = Ballp(N(3), N[4, 5], N(6))
    @test minkowski_sum(B1, B2) == Ballp(N(3), N[5, 7], N(9))
end
