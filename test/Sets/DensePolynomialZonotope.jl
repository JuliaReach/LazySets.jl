for N in [Float64, Float32, Rational{Int}]
    # example set from Figure 2 in original paper
    c = zeros(N, 2)
    E1 = Matrix(Diagonal(N[-1, 0.5]))
    E2 = N[1 1; 0.5 0.3]
    E = [E1, E2]
    F2 = hcat(N[-0.5, 1])
    F = [F2]
    G = hcat(N[0.3, 0.3])
    P = DensePolynomialZonotope(c, E, F, G)

    @test dim(P) == 2
    @test polynomial_order(P) == 2
    @test center(P) == c
    @test ngens_dep(P) == 5
    @test ngens_indep(P) == 1
    @test order(P) == 3 // 1

    # type-specific concrete methods
    # TODO these commands do not test anything
    P2 = copy(P)
    scale!(N(2), P2)
    @test P2 == scale(N(2), P)
    linear_map(N[1.0 2.0; 2.0 5.0], P)
    z = Zonotope(N[1.0, 2.0], Matrix(N(1)I, 2, 2))
    minkowski_sum(P, z)
    minkowski_sum(z, P)
end
