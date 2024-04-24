for N in [Float64, Float32, Rational{Int}]
    @test rand(SparsePolynomialZonotope) isa SparsePolynomialZonotope

    # Example 3.1.2 from thesis
    c = N[4, 4]
    G = N[2 1 2; 0 2 2]
    GI = hcat(N[1; 0])
    E = [1 0 3; 0 1 1]
    PZ = SparsePolynomialZonotope(c, G, GI, E)
    # Example 3.1.21 from thesis
    PZ2 = SparsePolynomialZonotope(zeros(N, 2), N[2 0 1; 1 2 1],
                                   zeros(N, 2, 0), [1 0 1; 0 1 3])

    @test center(PZ) == c
    @test genmat_dep(PZ) == G
    @test genmat_indep(PZ) == GI
    @test expmat(PZ) == E
    @test indexvector(PZ) == [1, 2]

    @test dim(PZ) == 2
    @test ngens_dep(PZ) == 3
    @test ngens_indep(PZ) == 1
    @test nparams(PZ) == 2
    @test order(PZ) == 2 // 1

    LM = linear_map(N[1//2 0; 0 1//2], PZ)
    @test center(LM) == [2, 2]
    @test genmat_dep(LM) == [1 0.5 1; 0 1 1]
    @test genmat_indep(LM) == hcat([0.5, 0.0])
    @test expmat(LM) == expmat(PZ)

    M = N[-0.5 0.2; -0.1 0.6]
    LMPZ = linear_map(M, PZ2)
    @test center(LMPZ) == zeros(N, 2)
    @test genmat_dep(LMPZ) ≈ N[-0.8 0.4 -0.3; 0.4 1.2 0.5] atol = 1e-7
    @test isempty(genmat_indep(LMPZ))
    @test expmat(LMPZ) == [1 0 1; 0 1 3]
    @test indexvector(LMPZ) == indexvector(PZ)

    ESPZ = PZ ⊞ PZ2
    @test center(ESPZ) == [4, 4]
    @test genmat_dep(ESPZ) == [2 1 2 2 0 1; 0 2 2 1 2 1]
    @test genmat_indep(ESPZ) == hcat([1, 0])
    @test expmat(ESPZ) == [1 0 3 1 0 1; 0 1 1 0 1 3]
    @test indexvector(ESPZ) == indexvector(PZ)

    TPZ = translate(PZ, N[1, 2])
    @test center(TPZ) == N[5, 6]
    @test genmat_dep(TPZ) == genmat_dep(TPZ)
    @test genmat_indep(TPZ) == genmat_indep(TPZ)
    @test expmat(TPZ) == expmat(TPZ)

    CPPZ = cartesian_product(PZ, PZ2)
    @test center(CPPZ) == N[4, 4, 0, 0]
    @test genmat_dep(CPPZ) == N[2 1 2 0 0 0;
                                0 2 2 0 0 0;
                                0 0 0 2 0 1;
                                0 0 0 1 2 1]
    @test genmat_indep(CPPZ) == hcat(N[1, 0, 0, 0])
    @test expmat(CPPZ) == [1 0 3 0 0 0;
                           0 1 1 0 0 0;
                           0 0 0 1 0 1;
                           0 0 0 0 1 3]
    Z = overapproximate(PZ2, Zonotope)
    CPPZ = cartesian_product(PZ, Z)
    @test center(CPPZ) == N[4, 4, 0, 0]
    @test genmat_dep(CPPZ) == N[2 1 2;
                                0 2 2;
                                0 0 0;
                                0 0 0]
    @test genmat_indep(CPPZ) == N[1 0 0 0;
                                  0 0 0 0;
                                  0 2 0 1;
                                  0 1 2 1]
    @test expmat(CPPZ) == [1 0 3;
                           0 1 1]

    MSPZ = minkowski_sum(PZ, PZ2)
    @test center(MSPZ) == [4, 4]
    @test genmat_dep(MSPZ) == [2 1 2 2 0 1; 0 2 2 1 2 1]
    @test genmat_indep(MSPZ) == hcat([1, 0])
    @test expmat(MSPZ) == [1 0 3 0 0 0;
                           0 1 1 0 0 0;
                           0 0 0 1 0 1;
                           0 0 0 0 1 3]

    S = SparsePolynomialZonotope(N[-0.5, -0.5], N[1.0 1 1 1; 1 0 -1 1], zeros(N, 2, 0),
                                 [1 0 1 2; 0 1 1 0])
    Z = overapproximate(S, Zonotope)
    @test isequivalent(Z, Zonotope(N[0, 0], N[1.5 1 1; 1.5 0 -1]))

    PZ = SparsePolynomialZonotope(N[-1, 2], N[1 2 0 2; 0 1 2 -1], N[1 0; 2 0], [1 0 1 2; 0 0 0 1])
    PZreduced = remove_redundant_generators(PZ)

    @test center(PZreduced) == N[1, 3]
    @test genmat_dep(PZreduced) == N[1 2; 2 -1]
    @test genmat_indep(PZreduced) == hcat(N[1, 2])
    @test expmat(PZreduced) == [1 2; 0 1]

    # support function (enclosure)
    for (d, v) in [(N[1, 0], N(3)), (N[1, 1], N(7)), (N[1, -1], N(3))]
        v1 = ρ(d, PZ2)  # default enclosure method
        v2 = ρ(d, PZ2; enclosure_method=RangeEnclosures.NaturalEnclosure())
        @test v <= v1 <= v2
    end
end

for N in [Float64]
    PZS = SimpleSparsePolynomialZonotope(N[0.2, -0.6], N[1 0; 0 0.4], [1 0; 0 1])
    PZ = convert(SparsePolynomialZonotope, PZS)
    @test center(PZ) == center(PZS)
    @test genmat_dep(PZ) == genmat(PZS)
    @test expmat(PZ) == expmat(PZS)
    @test isempty(genmat_indep(PZ))
    @test indexvector(PZ) == 1:2

    # conversion from Taylor model
    x₁, x₂, x₃ = set_variables(Float64, ["x₁", "x₂", "x₃"]; order=3)
    dom1 = IA.interval(N(-1), N(1))
    dom = dom1 × dom1 × dom1
    x0 = IA.IntervalBox(IA.mid.(dom)...)
    rem = IA.interval(N(0), N(0))
    p₁ = 33 + 2x₁ + 3x₂ + 4x₃ + 5x₁^2 + 6x₂ * x₃ + 7x₃^2 + 8x₁ * x₂ * x₃
    p₂ = x₃ - x₁
    vTM = [TaylorModels.TaylorModelN(pi, rem, x0, dom) for pi in [p₁, p₂]]
    PZ = convert(SparsePolynomialZonotope, vTM)
    # the following tests check for equality (but equivalence should be tested)
    @test PZ.c == N[33, 0]
    @test PZ.G == N[2 3 4 5 6 7 8;
                    -1 0 1 0 0 0 0]
    @test PZ.GI == Matrix{N}(undef, 2, 0)
    @test PZ.E == [1 0 0 2 0 0 1;
                   0 1 0 0 1 0 1;
                   0 0 1 0 1 2 1]
    # interestingly, the zonotope approximations are equivalent
    Zt = overapproximate(vTM, Zonotope)
    Zp = overapproximate(PZ, Zonotope)
    @test isequivalent(Zt, Zp)

    # conversion back to Taylor model
    vTM2 = convert(Vector{<:TaylorModelN}, PZ)
    @test vTM == vTM2
end

for Z in [rand(Zonotope), rand(Hyperrectangle)]
    ZS = convert(SparsePolynomialZonotope, Z)
    @test center(ZS) == center(Z)
    @test genmat_dep(ZS) == genmat(Z)
    @test isempty(genmat_indep(ZS))
    @test expmat(ZS) == I
end

let
    # example 3.1.32 from Niklas' thesis (page 60)
    PZ1 = SparsePolynomialZonotope([0.0, 0.0], [1.0 -1 1; -1 2 1], hcat([0.1, 0]), [1 0 2; 0 1 1],
                                   1:2)
    Q = QuadraticMap([[0.5 0.5; 1 -0.5], [-1.0 0; 1 0]], PZ1)
    QPZ = overapproximate(Q, SparsePolynomialZonotope)
    @test center(QPZ) ≈ [0.0025, -0.005]
    @test genmat_dep(QPZ) == [-1.5 5.5 2.0 -4.5 -1.5 1.5;
                              -2.0 5.0 -2.0 -3.0 3.0 0.0]
    @test genmat_indep(QPZ) ≈ [-0.05 0.2 0.25 0.0025;
                               -0.3 0.4 -0.10 -0.005]
    @test expmat(QPZ) == [2 1 3 0 2 4;
                          0 1 1 2 2 2]

    # example 3.1.40 from Nilkas' thesis (page 68)
    c = [0.0, 0]
    G = [-1.0 -2.0 -1.0 2.0 0.01 0.4
         1.0 0.0 -1.0 1.0 0.2 0.0]

    GI = [0.2 0.01
          0.02 -0.4]

    E = [1 0 1 2 2 0
         0 1 1 0 0 2
         0 0 0 0 1 2]

    P = SparsePolynomialZonotope(c, G, GI, E)
    Pred = reduce_order(P, 3)

    @test center(Pred) == [0.2, 0]
    @test genmat_dep(Pred) == [-1 -2 -1 2;
                               1 0 -1 1]

    @test genmat_indep(Pred) ≈ [0.42 0; 0 0.62]

    @test expmat(Pred) == [1 0 1 2; 0 1 1 0]
    @test indexvector(Pred) == [1, 2]
end
