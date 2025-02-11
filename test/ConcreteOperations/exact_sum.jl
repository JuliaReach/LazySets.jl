for N in [Float64, Float32, Rational{Int}]
    B = Ball1(N[0, 0], N(1))
    c = N[4, 4]
    G = N[2 1 2; 0 2 2]
    GI = hcat(N[1; 0])
    E = [1 0 3; 0 1 1]
    PZ = SparsePolynomialZonotope(c, G, GI, E)

    @test_throws ArgumentError exact_sum(B, PZ)
    @test_throws ArgumentError exact_sum(PZ, B)
end
