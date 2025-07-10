using LazySets, Test, LinearAlgebra
import IntervalArithmetic as IA
@static if VERSION >= v"1.9"
    vIA = pkgversion(IA)
else
    import PkgVersion
    vIA = PkgVersion.Version(IA)
end

for N in [Float64, Float32, Rational{Int}]
    # example from slide 13 of Niklas talk at JuliaReach & JuliaIntervals Days 3
    c = N[2, 0]
    G = N[1 2; 2 2.0]
    E = [1 4; 1 2]

    S = SimpleSparsePolynomialZonotope(c, G, E)

    @test genmat(S) == G
    @test expmat(S) == E
    @test dim(S) == 2
    @test ngens(S) == 2
    @test ngens_dep(S) == 2
    @test ngens_indep(S) == 0
    @test genmat_dep(S) == G
    @test genmat_indep(S) == zeros(N, 2, 0)
    @test nparams(S) == 2
    @test order(S) == 1 // 1
    @test polynomial_order(S) == 6

    Z = Zonotope(N[3.0, 1], N[1 1; 2 1.0])
    @test overapproximate(S, Zonotope) == Z
    @test overapproximate(S, UnionSetArray{Zonotope}; nsdiv=1) == UnionSetArray([Z])
    @test length(overapproximate(S, UnionSetArray{Zonotope}; nsdiv=3)) == 9
    if vIA >= v"0.19.0"  # `mince` with non-uniform partition was introduced in IA v0.19.0
        @test length(overapproximate(S, UnionSetArray{Zonotope}; partition=(2, 3))) == 6
    end

    LMS = linear_map(N[1 2; 3 4], S)
    @test LMS isa SimpleSparsePolynomialZonotope{N}
    @test center(LMS) == N[2, 6]
    @test genmat(LMS) == N[5 6; 11 14]
    @test expmat(LMS) == expmat(S)

    LMS2 = linear_map(N[1//2 0; 0 1//2], S)
    @test center(LMS2) == N[1, 0]
    @test genmat(LMS2) == N[0.5 1; 1 1]
    @test expmat(LMS2) == expmat(LMS2)

    MSS = minkowski_sum(S, S)
    @test center(MSS) == N[4, 0]
    @test genmat(MSS) == N[1 2 1 2; 2 2 2 2.0]
    @test expmat(MSS) == [1 4 0 0; 1 2 0 0; 0 0 1 4; 0 0 1 2]

    CPS = cartesian_product(S, S)
    @test center(CPS) == N[2, 0, 2, 0]
    @test genmat(CPS) == N[1 2 0 0; 2 2 0 0; 0 0 1 2; 0 0 2 2.0]
    @test expmat(CPS) == [1 4 0 0; 1 2 0 0; 0 0 1 4; 0 0 1 2]
    Z = overapproximate(S, Zonotope)
    PZZ = cartesian_product(S, Z)
    @test center(PZZ) == vcat(center(S), center(Z))
    @test genmat_dep(PZZ) == vcat(genmat_dep(S), zeros(N, 2, 2))
    @test genmat_indep(PZZ) == vcat(zeros(N, 2, 2), genmat(Z))
    @test expmat(PZZ) == expmat(S)

    _c = N[2, 0]
    _g = N[0 0.5 1 0.5 1 0.5 1 -0.5 -1
           0 1 1 1 1 1 1 -1 -1]

    _e = [0 1 4 1 4 0 0 0 0
          0 1 2 1 2 0 0 0 0
          0 0 0 0 0 1 4 1 4
          0 0 0 0 0 1 2 1 2
          1 0 0 1 1 0 0 1 1]

    LCS = linear_combination(S, S)

    @test center(LCS) == _c
    @test genmat(LCS) == _g
    @test expmat(LCS) == _e

    CH1 = convex_hull(S)

    @test center(CH1) == _c
    @test genmat(CH1) == _g
    @test expmat(CH1) == _e

    # convex hull with itself
    CH2 = convex_hull(S, S)
    Z1 = remove_redundant_generators(overapproximate(CH1, Zonotope))
    Z2 = remove_redundant_generators(overapproximate(CH2, Zonotope))
    @test_broken isequivalent(Z1, Z2)
    @test Z1 ⊆ Z2 && center(Z1) == center(Z2)

    # quadratic map
    S2 = SimpleSparsePolynomialZonotope(N[1 // 5, -3 // 5], N[1 0; 0 2//5], [1 0; 0 1])
    Q = [N[0 1; 1 0], N[-16//5 6//5; 6//5 0]]
    q1 = quadratic_map(Q, S2)
    @test concretize(QuadraticMap(Q, S2)) == q1
    q2 = quadratic_map(Q, S2, S2)

    @test center(q1) == center(q2) ≈ N[-6 // 25, -52 // 125]
    @test genmat(q1) ≈ genmat(q2) ≈ N[-6//5 4//25 0 4//5;
                                      -272//100 192//1000 -16//5 96//100]

    @test expmat(q1) == expmat(q2) == [1 0 2 1;
                                       0 1 0 1]

    c = N[1, 2, 3]

    G = N[0 1 4 7 10 13 16
          0 2 5 8 11 14 17
          0 3 6 9 12 15 18]

    E = [1 4 1 3 0 3 7
         2 5 2 4 0 4 8
         3 6 3 5 0 5 9]

    Sred = remove_redundant_generators(SimpleSparsePolynomialZonotope(c, G, E))
    @test center(Sred) == N[11, 13, 15]
    @test genmat(Sred) == N[1 4 20 16
                            2 5 22 17
                            3 6 24 18]
    @test expmat(Sred) == [4 1 3 7
                           5 2 4 8
                           6 3 5 9]

    # isoperationtype
    @test !isoperationtype(SimpleSparsePolynomialZonotope)

    # extrema approximation
    l1, u1 = extrema(S; algorithm="zonotope")
    @test (l1, u1) == extrema(S)  # default algorithm
    @test_throws ArgumentError extrema(S; algorithm="???")
    H1 = Hyperrectangle(; low=l1, high=u1)
    l2, u2 = extrema(S; algorithm="lowhigh")
    H2 = Hyperrectangle(; low=l2, high=u2)
    @test H2 ⊆ H1 == Hyperrectangle(N[3, 1], N[2, 3])

    # support function (enclosure)
    P = SimpleSparsePolynomialZonotope(zeros(N, 2), N[2 0 1; 1 2 1], [1 0 1; 0 1 3])
    for (d, v) in [(N[1, 0], N(3)), (N[1, 1], N(7)), (N[1, -1], N(3))]
        v1 = ρ(d, P)  # default enclosure method
        v2 = ρ(d, P; enclosure_method=RangeEnclosures.NaturalEnclosure())
        @test v <= v1 <= v2
    end

    # translate
    TPZ = translate(P, N[1, 2])
    @test center(TPZ) == N[1, 2]
    @test genmat_dep(TPZ) == genmat_dep(P)
    @test expmat(TPZ) == expmat(P)
end

for N in [Float64, Float32]
    # rand
    @test rand(SimpleSparsePolynomialZonotope; N=N) isa SimpleSparsePolynomialZonotope{N}
end

for Z in [rand(Zonotope), rand(Hyperrectangle)]
    ZS = convert(SimpleSparsePolynomialZonotope, Z)
    @test center(ZS) == center(Z)
    @test genmat(ZS) == genmat(Z)
    @test expmat(ZS) == I
end

let
    SPZ = SparsePolynomialZonotope([4.0, 4], [2.0 1 2; 0 2 2], hcat([1.0; 0]), [1 0 3; 0 1 1])
    SSPZ = convert(SimpleSparsePolynomialZonotope, SPZ)
    @test center(SSPZ) == center(SPZ)
    @test genmat(SSPZ) == hcat(SPZ.G, SPZ.GI)
    @test expmat(SSPZ) == [1 0 3 0; 0 1 1 0; 0 0 0 1]
end
