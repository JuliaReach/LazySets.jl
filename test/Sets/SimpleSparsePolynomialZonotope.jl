for N in [Float64, Float32, Rational{Int}]
    # example from slide 13 of Niklas talk at JuliaReach & JuliaIntervals Days 3
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

    Z = rand(Zonotope)
    ZS = convert(SimpleSparsePolynomialZonotope, Z)
    @test center(ZS) == center(Z)
    @test genmat(ZS) == genmat(Z)
    @test expmat(ZS) == I

    @test overapproximate(S, Zonotope) == Zonotope(N[3., 1], N[1 1;2 1.])
    @test length(overapproximate(S, Zonotope; nsdiv=3)) == 9
    @test length(overapproximate(S, Zonotope; partition=(2, 3))) == 6

    LMS = linear_map([1.0 2.0;3.0 4.0], S);
    @test center(LMS) == [2.0, 6.0]
    @test genmat(LMS) == [5.0 6.0;11.0 14.0]
    @test expmat(LMS) == expmat(S)

    MSS = minkowski_sum(S, S)
    @test center(MSS) == [4.0, 0.0]
    @test genmat(MSS) == [1 2 1 2;2 2 2 2.]
    @test expmat(MSS) == [1 4 0 0;1 2 0 0;0 0 1 4;0 0 1 2]

    CPS = cartesian_product(S, S)
    @test center(CPS) == [2, 0, 2.0, 0]
    @test genmat(CPS) == [1 2 0 0;2 2 0 0;0 0 1 2;0 0 2 2.]
    @test expmat(CPS) == [1 4 0 0;1 2 0 0;0 0 1 4;0 0 1 2]

    _c = N[2.0, 0.0]
    _g = N[ 0.0  0.5  1.0  0.5  1.0  0.5  1.0  -0.5  -1.0
           0.0  1.0  1.0  1.0  1.0  1.0  1.0  -1.0  -1.0]

    _e = [0  1  4  1  4  0  0  0  0
          0  1  2  1  2  0  0  0  0
          0  0  0  0  0  1  4  1  4
          0  0  0  0  0  1  2  1  2
          1  0  0  1  1  0  0  1  1]

    LCS = linear_combination(S, S)

    @test center(LCS) == _c
    @test genmat(LCS) == _g
    @test expmat(LCS) == _e

    CH1 = convex_hull(S)

    @test center(CH1) == _c
    @test genmat(CH1) == _g
    @test expmat(CH1) == _e
end
