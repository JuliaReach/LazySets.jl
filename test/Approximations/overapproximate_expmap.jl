
using Test, LazySets, LinearAlgebra
using IntervalMatrices: IntervalMatrix
using RangeEnclosures

for N in @tN([Float32, Float64, Rational{Int}])
    MZ = MatrixZonotope(N[1 1; -1 1], [N[1 0; 1 2]])
    P = SparsePolynomialZonotope(N[1, -1], N[1 1; 0 -1], hcat(N[0, 1]), [2 1; 0 1; 1 0], [1, 2, 3])

    # general test
    mzexp = MatrixZonotopeExp(MZ)
    em = ExponentialMap(mzexp, P)
    Pex = overapproximate(em, SparsePolynomialZonotope, 2)
    @test Pex isa SparsePolynomialZonotope{N}
    @test dim(Pex) == dim(P)
end

for N in @tN([Float32, Float64])
    P = SparsePolynomialZonotope(N[1, -1], N[1 1; 0 -1], hcat(N[0, 1]), [2 1; 0 1; 1 0], [1, 2, 3])

    # expectation: exp(0)*P = I*P = P
    MZ = MatrixZonotope(zeros(N, 2, 2), Vector{Matrix{N}}())
    mzexp = MatrixZonotopeExp(MZ)
    em = ExponentialMap(mzexp, P)
    Pex = overapproximate(em, SparsePolynomialZonotope, 2)
    @test center(Pex) == center(P)
    @test genmat_dep(Pex) == genmat_dep(P)
    @test genmat_indep(Pex) == genmat_indep(P)
    @test expmat(Pex) == expmat(P)

    MZ = MatrixZonotope(N[1 -1; -1 2], [N[1.001 -0.999; -0.999 2.005]])
    mzexp =MatrixZonotopeExp(MZ)
    em_spz = ExponentialMap(mzexp, P)
    Pex = overapproximate(em_spz,SparsePolynomialZonotope, 3)

    Z = overapproximate(convert(SimpleSparsePolynomialZonotope, P), Zonotope)
    em_z = ExponentialMap(mzexp, Z)
    Zex = overapproximate(em_z, Zonotope, 3)

    # samples from the SPZ and check if the points are contained in the overapproximation
    sampler = LazySets.PolynomialZonotopeSampler()
    pts = sample(Pex, 10; sampler=sampler)
    @test all(p ∈ Zex for p in pts)
end
