for N in [Float64, Float32, Rational{Int}]

    c = zeros(N, 2)
    E1 = Matrix(Diagonal(N[-1, 0.5]))
    E2 = N[1 1; 0.5 0.3]
    E = [E1, E2]
    F2 = hcat(N[-0.5, 1])
    F = [F2]
    G = Matrix(N(0.3)I, 2, 2)
    p = PolynomialZonotope(c, E, F, G)

    @test dim(p) == 2
    @test order(p) == 7//2
    @test polynomial_order(p) == 2

    # type-specific concrete methods
    scale(N(0.5), p)
    linear_map(N[1.0 2.0; 2.0 5.0], p)
    z = Zonotope(N[1., 2.], Matrix(N(1)I, 2, 2))
    minkowski_sum(p, z)
    minkowski_sum(z, p)
end