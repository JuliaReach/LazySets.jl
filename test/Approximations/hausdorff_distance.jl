_in_interval(v, x, ε) = x - ε <= v <= x + ε

for N in @tN([Float64, Float32, Rational{Int}])
    ε = N(1e-3)

    b1 = BallInf(zeros(N, 2), N(1))
    b2 = BallInf(ones(N, 2), N(1))

    dist = hausdorff_distance(b1, b2; p=N(Inf), ε=ε)
    @test _in_interval(dist, N(1), ε)

    dist = hausdorff_distance(b1, b1; p=N(Inf), ε=ε)
    @test _in_interval(dist, N(0), ε)
end

if test_suite_polyhedra
    for N in [Float64]
        ε = N(1e-3)

        b1 = N[1 0; 0 1] * BallInf(zeros(N, 2), N(1))
        b2 = BallInf(ones(N, 2), N(1))

        dist = hausdorff_distance(b1, b2; p=N(Inf), ε=ε)
        @test _in_interval(dist, N(1), ε)
    end
end
