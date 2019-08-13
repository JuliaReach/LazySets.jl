_in_interval(v, x, ε) = x - ε <= v <= x + ε

for N in [Float64, Float32, Rational{Int}]
    ε = N(1e-3)

    b1 = BallInf(zeros(N, 2), N(1))
    b1_lazy = N[1 0; 0 1] * b1
    b2 = BallInf(ones(N, 2), N(1))

    for b in [b1, b1_lazy]
        hd = hausdorff_distance(b, b2; p=N(Inf), ε=ε)
        @test _in_interval(hd, N(1), ε)
    end
end
