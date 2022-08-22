for N in [Float64, Float32, Rational{Int}]

    @test rand(SparsePolynomialZonotope) isa SparsePolynomialZonotope

    # example from Niklas thesis (page 32)
    c = N[4, 4]
    G = N[2 1 2;0 2 2]
    GI = hcat(N[1; 0])
    E = [1 0 3;0 1 1]
    PZ = SparsePolynomialZonotope(c, G, GI, E)

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


    LM = linear_map(0.5, PZ)
    @test center(LM) == [2, 2]
    @test genmat_dep(LM) == [1 0.5 1;0 1 1]
    @test genmat_indep(LM) == hcat([0.5, 0.0])
    @test expmat(LM) == expmat(PZ)

    M = N[-0.5 0.2;-0.1 0.6]
    PZ2 = SparsePolynomialZonotope(zeros(N, 2), N[2 0 1;1 2 1], zeros(N, 2, 0), [1 0 1;0 1 3])
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
    PZ1 = SparsePolynomialZonotope([0.0, 0.0], [1. -1 1;-1 2 1], hcat([0.1, 0]), [1 0 2;0 1 1], 1:2)
    Q = QuadraticMap([[0.5 0.5;1 -0.5], [-1.0 0;1 0]], PZ1)
    QPZ = overapproximate(Q, SparsePolynomialZonotope)
    @test_broken center(QPZ) == [0.0025, -0.005]
    @test_broken genmat_dep(QPZ) == [-4.5 5.5 -1.5 -1.5 2 1.5;-3 5 -2 3 -2 0]
    @test_broken genmat_indep(QPZ) == [0.0025 -0.05 0.15 0.15 0 0.05 0.1;-0.005 -0.2 0.3 0 -0.1 0.1 -0.1]
    @test_broken expmat(QPZ) == [0 1 2 2 3 4;2 1 0 2 1 2]
end
