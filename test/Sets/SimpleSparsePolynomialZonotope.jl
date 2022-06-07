for N in [Float64, Float32, Rational{Int}]
    # example from slide 13 of Niklas talk at JRJID3
    c = N[2.0, 0.0]
    G = N[1 2;2 2.]
    E = [1 4;1 2]

    S = SimpleSparsePolynomialZonotope(c, G, E)

    @test genmat(S) == G
    @test expmat(S) == E
    @test dim(S) == 2
    @test ngens(S) == 2
    @test nparams(S) == 2
    @test order(S) == 1 // 1

    @test overapproximate(S, Zonotope) == Zonotope(N[3., 1], N[1 1;2 1.])
    @test length(overapproximate(S, Zonotope; nsdiv=3)) == 9
    @test length(overapproximate(S, Zonotope; partition=(2, 3))) == 6
end
