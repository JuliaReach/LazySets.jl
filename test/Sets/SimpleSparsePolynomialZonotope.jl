for N in [Float64, Float32, Rational{Int}]
    # example from slide 13 of Niklas talk at JuliaReach & JuliaIntervals Days 3
    c = N[2.0, 0.0]
    G = N[1 2;2 2.]
    E = [1 4;1 2]

    @test rand(SimpleSparsePolynomialZonotope) isa SimpleSparsePolynomialZonotope

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

    LMS = linear_map([1.0 2.0;3.0 4.0], S);
    @test center(LMS) == [2.0, 6.0]
    @test genmat(LMS) == [5.0 6.0;11.0 14.0]
    @test expmat(LMS) == expmat(S)

    LMS2 = linear_map(0.5, S)
    @test center(LMS2) == [1.0, 0.0]
    @test genmat(LMS2) == [0.5 1;1 1]
    @test expmat(LMS2) == expmat(LMS2)

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

    S2 = SimpleSparsePolynomialZonotope(N[0.2, -0.6], N[1 0;0 0.4], [1 0;0 1])
    Q = [[0. 1;1 0], [-3.2 1.2;1.2 0]]

    q1 = quadratic_map(Q, S2)

    @test center(q1) ≈ [-0.24, -0.416]
    @test genmat(q1) ≈ [-1.2   0.16    0.0  0.4   0.4   0.0;
                         -2.72  0.192  -3.2  0.48  0.48  0.0]

    @test expmat(q1) == [ 1  0  2  1  1  0;
                          0  1  0  1  1  2]

    c = N[1., 2, 3]

    G = N[0.0  1.0  4.0  7.0  10.0  13.0  16.0
         0.0  2.0  5.0  8.0  11.0  14.0  17.0
         0.0  3.0  6.0  9.0  12.0  15.0  18.0]

    E = [1  4  1  3  0  3  7
         2  5  2  4  0  4  8
         3  6  3  5  0  5  9]

    Sred = remove_redundant_generators(SimpleSparsePolynomialZonotope(c, G, E))
    @test center(Sred) == [11., 13, 15]
    @test genmat(Sred) == N[ 1.0  4.0  20.0  16.0
                            2.0  5.0  22.0  17.0
                            3.0  6.0  24.0  18.0]

    @test expmat(Sred) == [ 4  1  3  7
                            5  2  4  8
                            6  3  5  9]
end


for Z in [rand(Zonotope), rand(Hyperrectangle)]
    ZS = convert(SimpleSparsePolynomialZonotope, Z)
    @test center(ZS) == center(Z)
    @test genmat(ZS) == genmat(Z)
    @test expmat(ZS) == I
end

SPZ = SparsePolynomialZonotope([4., 4], [2. 1 2;0 2 2], hcat([1.; 0]), [1 0 3;0 1 1])
SSPZ = convert(SimpleSparsePolynomialZonotope, SPZ)
@test center(SSPZ) == center(SPZ)
@test genmat(SSPZ) == hcat(SPZ.G, SPZ.GI)
@test expmat(SSPZ) == [1 0 3 0;0 1 1 0;0 0 0 1]
