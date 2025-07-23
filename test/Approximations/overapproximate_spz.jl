using Test, LazySets, LinearAlgebra
using IntervalMatrices: IntervalMatrix
using RangeEnclosures
using IntervalArithmetic: Interval

for N in @tN([Float32, Float64, Rational{Int}])
    # test linear map
    P = SparsePolynomialZonotope(N[1, -1], N[1 1; 0 -1], hcat(N[0, 1]), [2 1; 0 1; 1 0], [1, 2, 3])
    MZ = MatrixZonotope(N[1 1; -1 1], [N[1 0; 1 2]])

    res = overapproximate(MZ * P, SparsePolynomialZonotope)
    @test center(res) == [0, -2]
    @test genmat_dep(res) == hcat(N[1 0; -1 -2], N[1, -1], N[1 1; 1 -1])
    @test genmat_indep(res) == hcat([1, 1], [0, 2])
    @test expmat(res) == hcat([2 1; 0 1; 1 0], [1; 0; 0], [3 2; 0 1; 1 0])

    # case: 0 gens matrix zonotope
    MZ = MatrixZonotope(N[1 1; -1 1], Vector{Matrix{N}}())
    res = overapproximate(MZ * P, SparsePolynomialZonotope)
    @test res == linear_map(N[1 1; -1 1], P)

    #case: id_mz ≠ id_spz
    MZ = MatrixZonotope(N[1 1; -1 1], [N[1 0; 1 2]], [4])
    res = overapproximate(MZ * P, SparsePolynomialZonotope)
    @test center(res) == [0, -2]
    @test genmat_dep(res) == hcat(N[1 0; -1 -2], N[1, -1], N[1 1; 1 -1])
    @test genmat_indep(res) == hcat([1, 1], [0, 2])
    @test expmat(res) == hcat([2 1; 0 1; 1 0; 0 0], [0, 0, 0, 1], [2 1; 0 1; 1 0; 1 1])

    #case: dimension mismatching
    MZ_err = rand(MatrixZonotope; dim=(3, 3))
    @test_throws AssertionError res = overapproximate(MZ_err * P, SparsePolynomialZonotope)

    #case: matrix zonotope product
    MZ2 = MatrixZonotope(N[1.1 0.9; -1.1 1.1], [N[1.1 -0.1; 0.9 2.1]])
    res = overapproximate(MZ*MZ2*P, SparsePolynomialZonotope)  
    P_in = overapproximate(MZ2*P, SparsePolynomialZonotope)
    @test res == overapproximate( MZ*P_in, SparsePolynomialZonotope)

    #exponential map
    em = ExponentialMap(MZ, P) #test MZP
    Pex = overapproximate(em, SparsePolynomialZonotope, 2)
    @test isa(Pex, SparsePolynomialZonotope)
    @test size(center(Pex), 1) == size(center(P), 1)

    # expectation: exp(0)*P = I*P = P
    ZM = MatrixZonotope(zeros(N, 2, 2), Vector{Matrix{N}}())
    em = ExponentialMap(ZM, P)
    Pex = overapproximate(em, SparsePolynomialZonotope, 2)
    @test center(Pex) == center(P)
    @test genmat_dep(Pex) == genmat_dep(P)
    @test genmat_indep(Pex) == genmat_indep(P)
    @test expmat(Pex) == expmat(P)
end

for N in @tN([Float32, Float64])
    P = SparsePolynomialZonotope(N[1, -1], N[1 1; 0 -1], hcat(N[0, 1]), [2 1; 0 1; 1 0], [1, 2, 3])
    
    MZ = MatrixZonotope(N[1 -1; -1 2], [N[1.001 -0.999; -0.999 2.005]]) #SPD 
    em_spz = ExponentialMap(MZ, P)
    Pex = overapproximate(em_spz,SparsePolynomialZonotope, 3)

    Z = overapproximate(convert(SimpleSparsePolynomialZonotope, P), Zonotope)
    em_z = ExponentialMap(MZ, Z)
    Zex = overapproximate(em_z, Zonotope, 3)

    # samples from the SPZ and check if the points are contained in the overapproximation
    sampler = LazySets.PolynomialZonotopeSampler()
    pts = sample(Pex, 100; sampler=sampler)
    @test all(p ∈ Zex for p in pts)
end