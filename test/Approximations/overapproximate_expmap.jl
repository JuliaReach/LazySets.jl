using LazySets, Test, SparseArrays
using LazySets.MatrixZonotopeModule: vectorize
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

for N in @tN([Float64, Float32, Rational{Int}])
    MZ = MatrixZonotope(N[1 1; -1 1], [N[1 0; 1 2]])
    P = SparsePolynomialZonotope(N[1, -1], N[1 1; 0 -1], hcat(N[0, 1]), [2 1; 0 1; 1 0], [1, 2, 3])

    @static if isdefined(@__MODULE__, :IntervalMatrices)
        # general test
        mzexp = MatrixZonotopeExp(MZ)
        em = ExponentialMap(mzexp, P)
        Pex = overapproximate(em, SparsePolynomialZonotope, 2)
        @test Pex isa SparsePolynomialZonotope{N}
        @test dim(Pex) == dim(P)
    end
end

for N in @tN([Float32, Float64])
    P = SparsePolynomialZonotope(N[1, -1], N[1 1; 0 -1], hcat(N[0, 1]), [2 1; 0 1; 1 0], [1, 3, 5])
    M = N[1 1; -1 1]

    @static if isdefined(@__MODULE__, :IntervalMatrices)
        # expectation: exp(0)*P = I*P = P
        MZ = MatrixZonotope(zeros(N, 2, 2), Vector{Matrix{N}}())
        mzexp = MatrixZonotopeExp(MZ)
        em = ExponentialMap(mzexp, P)
        Pex = overapproximate(em, SparsePolynomialZonotope, 2)
        @test Pex == P

        MZ = MatrixZonotope(M, [N[1.001 -0.999; -0.999 2.005]])
        mzexp = MatrixZonotopeExp(MZ)
        em_spz = ExponentialMap(mzexp, P)
        Pex = overapproximate(em_spz, SparsePolynomialZonotope, 5)

        matnorm = norm(MZ, Inf) #in this case the exact norm is the same as the approximated up to rounding errors
        Pexnorm = overapproximate(em_spz, SparsePolynomialZonotope, 5; matnorm=matnorm)
        @test isapprox(center(Pexnorm), center(Pex); atol=1e-2)
        @test isapprox(genmat_dep(Pexnorm), genmat_dep(Pex); atol=1e-2)
        @test isapprox(genmat_indep(Pexnorm), genmat_indep(Pex); atol=1e-2)

        Z = overapproximate(convert(SimpleSparsePolynomialZonotope, P), Zonotope)
        em_z = ExponentialMap(mzexp, Z)
        Zex = overapproximate(em_z, Zonotope, 5)

        # samples from the SPZ and check if the points are contained in the overapproximation
        sampler = LazySets.PolynomialZonotopeSampler()
        pts = sample(Pex, 10; sampler=sampler)
        @test all(p ∈ Zex for p in pts)

        # overapproximate MatrixZonotopeExp
        
        # test degenerate case
        c = N[1 -2; 2 -1]
        A = MatrixZonotope(c, [zeros(N, 2, 2)])
        expA = MatrixZonotopeExp(A)
        res = overapproximate(expA, MatrixZonotope, 20) #for large k it converges to exp(c)
        @test isapprox(center(res), exp(c))
        @test ngens(res)== 0

        # degenerate case + product 
        B = MatrixZonotope(N[1 0; 0 1], Matrix{N}[])
        expAB = MatrixZonotopeExp(A * B)
        res2 = overapproximate(expA, MatrixZonotope, 20)
        @test isapprox(center(res), center(res2))

        @static if isdefined(@__MODULE__, :ExponentialUtilities) || isdefined(@__MODULE__, :Expokit)
            # matrix
            mzexp = SparseMatrixExp(sparse(M))
            em = ExponentialMap(mzexp, P)
            MPex = overapproximate(em, SparsePolynomialZonotope, 2)
            MZex = overapproximate(MPex, Zonotope)
            @test MZex ⊆ Zex
        end
    end
end

for N in @tN([Float64])
    # This test is only run with Float64 since it is expensive
    # and requires high Taylor order `k` to pass on Float32
    @static if isdefined(@__MODULE__, :IntervalMatrices)
        #inclusion 
        C = MatrixZonotope(N[1 -2; 2 -1], [N[0.1 0.05; 0 0.1]])
        expC = MatrixZonotopeExp(C)
        res = reduce_order(overapproximate(expC, MatrixZonotope, 4), 1)
        res2 = reduce_order(overapproximate(expC, MatrixZonotope, 8), 1)
        @test vectorize(res2) ⊆ vectorize(res)       
    end
end