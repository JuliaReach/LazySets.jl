for N in [Float64, Float32, Rational{Int}]
    X = Interval(N(1), N(2))
    Y = X + X
    if test_suite_polyhedra
        minkowski_sum(X, Y)
        minkowski_sum(X, Y; algorithm=Polyhedra.FourierMotzkin())
    end
end
