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
end

for N in [Float64, Float32]
    B1 = Ball2(N[1, 2], N(3))
    B2 = Ball2(N[4, 5], N(6))
    @test minkowski_sum(B1, B2) == Ball2(N[5, 7], N(9))

    B1 = Ballp(N(3), N[1, 2], N(3))
    B2 = Ballp(N(3), N[4, 5], N(6))
    @test minkowski_sum(B1, B2) == Ballp(N(3), N[5, 7], N(9))
end
