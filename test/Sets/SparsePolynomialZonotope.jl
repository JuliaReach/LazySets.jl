using LazySets, Test, LinearAlgebra
import IntervalArithmetic as IA
using IntervalArithmetic: IntervalBox
using LazySets.SparsePolynomialZonotopeModule: merge_id
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

for N in @tN([Float64, Float32, Rational{Int}])
    # Example 3.1.2 from thesis
    c = N[4, 4]
    G = N[2 1 2; 0 2 2]
    GI = hcat(N[1; 0])
    E = [1 0 3; 0 1 1]
    PZ = SparsePolynomialZonotope(c, G, GI, E)
    # Example 3.1.21 from thesis
    PZ2 = SparsePolynomialZonotope(zeros(N, 2), N[2 0 1; 1 2 1], zeros(N, 2, 0), [1 0 1; 0 1 3])

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
    @test LM isa SparsePolynomialZonotope{N}
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

    # exact sum: same IDs
    ESPZ = exact_sum(PZ, PZ2)
    @test center(ESPZ) == [4, 4]
    @test genmat_dep(ESPZ) == [2 1 2 2 0 1; 0 2 2 1 2 1]
    @test genmat_indep(ESPZ) == hcat([1, 0])
    @test expmat(ESPZ) == [1 0 3 1 0 1; 0 1 1 0 1 3]
    @test indexvector(ESPZ) == indexvector(PZ)
    #exact sum: different IDs
    PZ3 = SparsePolynomialZonotope(N[1, -1], N[1 -1; 0 2],
                                   hcat(N[0; 1]), [1 0; 2 1], [2, 3])
    ESPZ2 = exact_sum(PZ, PZ3)
    @test center(ESPZ2) == N[5, 3]
    @test genmat_dep(ESPZ2) == N[2 1 2 1 -1; 0 2 2 0 2]
    @test genmat_indep(ESPZ2) == N[1 0; 0 1]
    @test expmat(ESPZ2) == [1 0 3 0 0; 0 1 1 1 0; 0 0 0 2 1]
    @test indexvector(ESPZ2) == [1, 2, 3]

    # translate / translate!
    TPZ = translate(PZ, N[1, 2])
    @test TPZ isa SparsePolynomialZonotope{N}
    @test center(TPZ) == N[5, 6]
    @test genmat_dep(TPZ) == genmat_dep(PZ)
    @test genmat_indep(TPZ) == genmat_indep(PZ)
    @test expmat(TPZ) == expmat(PZ)
    @test indexvector(TPZ) == indexvector(PZ)
    P2 = copy(PZ)
    translate!(P2, N[1, 2])
    @test P2 == TPZ

    # cartesian_product SPZ/SPZ
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
    # cartesian_product mixed: SPZ/SSPZ
    CPPZ2 = cartesian_product(PZ, convert(SimpleSparsePolynomialZonotope, PZ2))
    @test CPPZ == CPPZ2
    CPPZ2 = cartesian_product(convert(SimpleSparsePolynomialZonotope, PZ), PZ2)
    @test center(CPPZ) == center(CPPZ2)  # the sets are equivalent, but this is hard to check
    # cartesian_product SPZ/Z
    Z = overapproximate(PZ2, Zonotope)
    SSPZ2 = convert(SimpleSparsePolynomialZonotope, PZ2)
    dom2 = IntervalBox(IA.interval(N(-1), N(1)), IA.interval(N(-1), N(1)))
    Z2 = overapproximate(SSPZ2, Zonotope, dom2)
    @test isequivalent(Z, Z2)
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

    # minkowski_sum SPZ/SPZ
    MSPZ = minkowski_sum(PZ, PZ2)
    @test center(MSPZ) == [4, 4]
    @test genmat_dep(MSPZ) == [2 1 2 2 0 1; 0 2 2 1 2 1]
    @test genmat_indep(MSPZ) == hcat([1, 0])
    @test expmat(MSPZ) == [1 0 3 0 0 0;
                           0 1 1 0 0 0;
                           0 0 0 1 0 1;
                           0 0 0 0 1 3]
    # minkowski_sum mixed SPZ/SSPZ
    MSPZ2 = minkowski_sum(PZ, convert(SimpleSparsePolynomialZonotope, PZ2))
    @test MSPZ == MSPZ2
    MSPZ2 = minkowski_sum(convert(SimpleSparsePolynomialZonotope, PZ), PZ2)
    @test center(MSPZ) == center(MSPZ2)  # the sets are equivalent, but this is hard to check

    S = SparsePolynomialZonotope(N[-0.5, -0.5], N[1.0 1 1 1; 1 0 -1 1], zeros(N, 2, 0),
                                 [1 0 1 2; 0 1 1 0])
    Z = overapproximate(S, Zonotope)
    @test isequivalent(Z, Zonotope(N[0, 0], N[1.5 1 1; 1.5 0 -1]))

    α = 2.0
    PZscaled = scale(α, PZ)
    @test center(PZscaled) == N[8, 8]
    @test genmat_dep(PZscaled) == N[4 2 4; 0 4 4]
    @test genmat_indep(PZscaled) == hcat(N[2, 0])
    @test expmat(PZscaled) == expmat(PZ)
    @test PZscaled == linear_map(α * Matrix{N}(I, 2, 2), PZ)

    # remove_redundant_generators
    PZ = SparsePolynomialZonotope(N[-1, 2], N[1 2 0 2; 0 1 2 -1], N[1 0; 2 0], [1 0 1 2; 0 0 0 1])
    PZreduced = remove_redundant_generators(PZ)
    @test center(PZreduced) == N[1, 3]
    @test genmat_dep(PZreduced) == N[1 2; 2 -1]
    @test genmat_indep(PZreduced) == hcat(N[1, 2])
    @test expmat(PZreduced) == [1 2; 0 1]
    # also removes almost-zero columns
    if N <: AbstractFloat
        PZ = SparsePolynomialZonotope(N[-1, 2], hcat(N[1e-16, -1e-16]),
                                      hcat(N[1e-16, -1e-16]), hcat([1, 1]))
        PZreduced = remove_redundant_generators(PZ)
        @test center(PZreduced) == N[-1, 2]
        @test size(genmat_dep(PZreduced)) == (2, 0)
        @test size(genmat_indep(PZreduced)) == (2, 0)
        @test size(expmat(PZreduced)) == (2, 0)
    end

    @static if isdefined(@__MODULE__, :RangeEnclosures)
        # support function (enclosure)
        for (d, v) in [(N[1, 0], N(3)), (N[1, 1], N(7)), (N[1, -1], N(3))]
            v1 = ρ(d, PZ2)  # default enclosure method
            v2 = ρ(d, PZ2; enclosure_method=RangeEnclosures.NaturalEnclosure())
            @test v <= v1 <= v2
        end

        # extrema approximation
        PZ = SparsePolynomialZonotope(N[-1, 2], N[1 2 0 2; 0 1 2 -1], N[1 0; 2 0],
                                      [1 0 1 2; 0 0 0 1])
        l1, u1 = extrema(PZ; algorithm="zonotope")
        @test (l1, u1) == extrema(PZ)  # default algorithm
        @test_throws ArgumentError extrema(PZ; algorithm="???")
        H1 = Hyperrectangle(; low=l1, high=u1)
        l2, u2 = extrema(PZ; algorithm="lowhigh")
        H2 = Hyperrectangle(; low=l2, high=u2)
        @test H2 ⊆ H1 == Hyperrectangle(N[0, 5 // 2], N[5, 11 // 2])
        # another example
        PZ = SparsePolynomialZonotope(N[-1 / 2, -1 / 2], N[1 1 1 1; 1 0 -1 1],
                                      Matrix{N}(undef, 2, 0), [1 0 1 2; 0 1 1 0])
        l1, u1 = extrema(PZ; algorithm="zonotope")
        H1 = Hyperrectangle(; low=l1, high=u1)
        l2, u2 = extrema(PZ; algorithm="lowhigh")
        H2 = Hyperrectangle(; low=l2, high=u2)
        @test H2 ⊆ H1 == Hyperrectangle(N[0, 0], N[7 // 2, 5 // 2])
        Z = Zonotope(N[0, 0], N[3//2 1 1; 3//2 0 -1])
        @test isequivalent(overapproximate(PZ, Zonotope), Z)
        SSPZ2 = convert(SimpleSparsePolynomialZonotope, PZ)
        @test isequivalent(overapproximate(SSPZ2, Zonotope, dom2), Z)
    end

    # linear_combination
    # example from slide 13 of Niklas talk at JuliaReach & JuliaIntervals Days 3
    c = N[2, 0]
    G = N[1 2; 2 2.0]
    E = [1 4; 1 2]
    SP = SimpleSparsePolynomialZonotope(c, G, E)
    P = convert(SparsePolynomialZonotope, SP)
    for P2 in (linear_combination(P, P), linear_combination(P, SP), linear_combination(SP, P))
        @test P2 isa SimpleSparsePolynomialZonotope{N}
        @test center(P2) == N[2, 0]
        @test genmat(P2) == N[0 0.5 1 0.5 1 0.5 1 -0.5 -1
                              0 1 1 1 1 1 1 -1 -1]
        @test expmat(P2) == [0 1 4 1 4 0 0 0 0
                             0 1 2 1 2 0 0 0 0
                             0 0 0 0 0 1 4 1 4
                             0 0 0 0 0 1 2 1 2
                             1 0 0 1 1 0 0 1 1]
    end
end

for N in @tN([Float64, Float32])
    # rand
    @test rand(SparsePolynomialZonotope; N=N) isa SparsePolynomialZonotope{N}
end

for Z in [rand(Zonotope), rand(Hyperrectangle)]
    ZS = convert(SparsePolynomialZonotope, Z)
    @test center(ZS) == center(Z)
    @test genmat_dep(ZS) == genmat(Z)
    @test isempty(genmat_indep(ZS))
    @test expmat(ZS) == I
end

let
    @static if isdefined(@__MODULE__, :TaylorModels)
        N = Float64
        PZS = SimpleSparsePolynomialZonotope(N[0.2, -0.6], N[1 0; 0 0.4], [1 0; 0 1])
        PZ = convert(SparsePolynomialZonotope, PZS)
        @test center(PZ) == center(PZS)
        @test genmat_dep(PZ) == genmat(PZS)
        @test expmat(PZ) == expmat(PZS)
        @test isempty(genmat_indep(PZ))
        @test indexvector(PZ) == 1:2

        # conversion from Taylor model
        x₁, x₂, x₃ = TaylorModels.set_variables(Float64, ["x₁", "x₂", "x₃"]; order=3)
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
        vTM2 = convert(Vector{<:TaylorModels.TaylorModelN}, PZ)
        @test vTM == vTM2
    end

    # isoperationtype
    @test !isoperationtype(SparsePolynomialZonotope)

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

    @test order(P) == 4
    for r in 4:8
        @test reduce_order(P, r) == P
    end

    # merge_id
    #case 0: dimensions mismatch
    id1 = [1, 2, 3]
    id2 = [1, 2]
    E = rand(Int, 2, 2)
    @test_throws ArgumentError merge_id(id1, id2, E, E)

    #case 1: identical IDs
    id = [10, 20, 30]
    E1 = rand(Int, 3, 4)
    E2 = rand(Int, 3, 2)

    Ē₁, Ē₂, idx = merge_id(id, id, E1, E2)

    @test idx == id
    @test Ē₁ === E1
    @test Ē₂ === E2

    # case 2: IDs with overlap
    E1 = [1 2; 1 0]
    id1 = [1, 2]

    E2 = [1 0 1; 3 2 0]
    id2 = [2, 3]

    Ē₁, Ē₂, idx = merge_id(id1, id2, E1, E2)

    @test idx == [1, 2, 3]
    @test Ē₁ == [1 2; 1 0; 0 0]
    @test Ē₂ == [0 0 0; 1 0 1; 3 2 0]

    #case 3: different IDs
    id2 = [3, 4]

    Ē₁, Ē₂, idx = merge_id(id1, id2, E1, E2)

    @test idx == [1, 2, 3, 4]
    @test Ē₁ == [1 2; 1 0; 0 0; 0 0]
    @test Ē₂ == [0 0 0; 0 0 0; 1 0 1; 3 2 0]
end
