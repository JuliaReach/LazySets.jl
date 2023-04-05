for N in [Float64, Float32, Rational{Int}]

    @test rand(SparsePolynomialZonotope) isa SparsePolynomialZonotope

    # Example 3.1.2 from thesis
    c = N[4, 4]
    G = N[2 1 2;0 2 2]
    GI = hcat(N[1; 0])
    E = [1 0 3;0 1 1]
    PZ = SparsePolynomialZonotope(c, G, GI, E)
    # Example 3.1.21 from thesis
    PZ2 = SparsePolynomialZonotope(zeros(N, 2), N[2 0 1;1 2 1],
                                   zeros(N, 2, 0), [1 0 1;0 1 3])

    @test center(PZ) == c
    @test genmat_dep(PZ) == G
    @test genmat_indep(PZ) == GI
    @test expmat(PZ) == E
    @test indexvector(PZ) == [1, 2]

    @test dim(PZ) == 2
    @test ngens_dep(PZ) == 3
    @test ngens_indep(PZ) == 1
    @test nparams(PZ) == 2
    @test order(PZ) == 2//1


    LM = linear_map(N[1//2 0; 0 1//2], PZ)
    @test center(LM) == [2, 2]
    @test genmat_dep(LM) == [1 0.5 1;0 1 1]
    @test genmat_indep(LM) == hcat([0.5, 0.0])
    @test expmat(LM) == expmat(PZ)

    M = N[-0.5 0.2;-0.1 0.6]
    LMPZ = linear_map(M, PZ2)
    @test center(LMPZ) == zeros(N, 2)
    @test genmat_dep(LMPZ) ≈ N[-0.8 0.4 -0.3;0.4 1.2 0.5] atol=1e-7
    @test isempty(genmat_indep(LMPZ))
    @test expmat(LMPZ) == [1 0 1;0 1 3]
    @test indexvector(LMPZ) == indexvector(PZ)

    ESPZ = PZ ⊞ PZ2
    @test center(ESPZ) == [4, 4]
    @test genmat_dep(ESPZ) == [2 1 2 2 0 1;0 2 2 1 2 1]
    @test genmat_indep(ESPZ) == hcat([1, 0])
    @test expmat(ESPZ) == [1 0 3 1 0 1;0 1 1 0 1 3]
    @test indexvector(ESPZ) == indexvector(PZ)

    TPZ = translate(PZ, N[1, 2])
    @test center(TPZ) == N[5, 6]
    @test genmat_dep(TPZ) == genmat_dep(TPZ)
    @test genmat_indep(TPZ) == genmat_indep(TPZ)
    @test expmat(TPZ) == expmat(TPZ)

    MSPZ = minkowski_sum(PZ, PZ2)
    @test center(MSPZ) == [4, 4]
    @test genmat_dep(MSPZ) == [2 1 2 2 0 1;0 2 2 1 2 1]
    @test genmat_indep(MSPZ) == hcat([1, 0])
    @test expmat(MSPZ) == [ 1  0  3  0  0  0;
                            0  1  1  0  0  0;
                            0  0  0  1  0  1;
                            0  0  0  0  1  3]

    S = SparsePolynomialZonotope(N[-0.5, -0.5], N[1. 1 1 1;1 0 -1 1], zeros(N, 2, 0), [1 0 1 2;0 1 1 0])
    Z = overapproximate(S, Zonotope)

    @test center(Z) == N[0, 0]
    @test genmat(Z) == N[1 1 1 0.5;1 0 -1 0.5]

    PZ = SparsePolynomialZonotope(N[-1, 2], N[1 2 0 2;0 1 2 -1], N[1 0;2 0], [1 0 1 2;0 0 0 1])
    PZreduced = remove_redundant_generators(PZ)

    @test center(PZreduced) == N[1, 3]
    @test genmat_dep(PZreduced) == N[1 2;2 -1]
    @test genmat_indep(PZreduced) == hcat(N[1, 2])
    @test expmat(PZreduced) == [1 2;0 1]

    # support function (enclosure)
    for (d, v) in [(N[1, 0], N(3)), (N[1, 1], N(7)), (N[1, -1], N(3))]
        v1 = ρ(d, PZ2)  # default enclosure method
        v2 = ρ(d, PZ2; enclosure_method=RangeEnclosures.NaturalEnclosure())
        @test v <= v1 <= v2
    end
end

SSPZ = SimpleSparsePolynomialZonotope([0.2, -0.6], [1 0;0 0.4], [1 0;0 1])
SPZ = convert(SparsePolynomialZonotope, SSPZ)
@test center(SPZ) == center(SSPZ)
@test genmat_dep(SPZ) == genmat(SSPZ)
@test expmat(SPZ) == expmat(SSPZ)
@test isempty(genmat_indep(SPZ))
@test indexvector(SPZ) == 1:2

for Z in [rand(Zonotope), rand(Hyperrectangle)]
    ZS = convert(SparsePolynomialZonotope, Z)
    @test center(ZS) == center(Z)
    @test genmat_dep(ZS) == genmat(Z)
    @test isempty(genmat_indep(ZS))
    @test expmat(ZS) == I
end

let
    # example 3.1.32 from Niklas' thesis (page 60)
    PZ1 = SparsePolynomialZonotope([0.0, 0.0], [1. -1 1;-1 2 1], hcat([0.1, 0]), [1 0 2;0 1 1], 1:2)
    Q = QuadraticMap([[0.5 0.5;1 -0.5], [-1.0 0;1 0]], PZ1)
    QPZ = overapproximate(Q, SparsePolynomialZonotope)
    @test center(QPZ) ≈ [0.0025, -0.005]
    @test genmat_dep(QPZ) ==  [-1.5  5.5   2.0  -4.5  -1.5  1.5;
                               -2.0  5.0  -2.0  -3.0   3.0  0.0]
    @test genmat_indep(QPZ) ≈ [-0.05 0.2  0.25  0.0025;
                               -0.3  0.4 -0.10 -0.005]
    @test expmat(QPZ) == [2  1  3  0  2  4;
                          0  1  1  2  2  2]


    # example 3.1.40 from Nilkas' thesis (page 68)
    c = [0.0, 0]
    G = [ -1.0  -2.0  -1.0  2.0  0.01  0.4
           1.0   0.0  -1.0  1.0  0.2   0.0]

    GI = [ 0.2    0.01
           0.02  -0.4]

    E = [ 1  0  1  2  2  0
          0  1  1  0  0  2
          0  0  0  0  1  2]

    P = SparsePolynomialZonotope(c, G, GI, E)
    Pred = reduce_order(P, 3)

    @test center(Pred) == [0.2, 0]
    @test genmat_dep(Pred) == [-1 -2 -1 2;
                                       1  0 -1 1]

    @test genmat_indep(Pred) ≈ [0.42 0;0 0.62]

    @test expmat(Pred) == [1 0 1 2;0 1 1 0]
    @test indexvector(Pred) == [1, 2]
end
