using LazySets, Test, SparseArrays
using LazySets.MatrixZonotopeModule: vectorize

if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

for N in @tN([Float32, Float64, Rational{Int}])
    #overapproximate matrix zonotope multiplication
    A = MatrixZonotope(N[1 1; -1 1], [N[1 0; 1 2]])
    B = MatrixZonotope(N[2 0; 1 -1], [N[0 1; 1 -1]])
    res = overapproximate(A * B, MatrixZonotope)

    ex_res = MatrixZonotope(N[3 -1; -1 -1], [N[0 1; 2 -1], N[2 0; 4 -2], N[1 0; 1 -2]])
    @static if isdefined(@__MODULE__, :Polyhedra)
        @test isequivalent(vectorize(res), vectorize(ex_res))
    end

    C = MatrixZonotope(N[1 0; 0 -1], [N[0 0; 1 0]])
    res2 = overapproximate(A * B * C, MatrixZonotope)
    @test res2 isa MatrixZonotope

    #empty generator
    D = MatrixZonotope(N[1 0; 0 -1], Matrix{N}[])
    res3 = overapproximate(A * D, MatrixZonotope)
    @test center(res3) == N[1 -1; -1 -1]
    @test generators(res3) == [N[1 0; 1 -2]]

    # more generators
    A = MatrixZonotope(N[1 2; 0 -1], [N[1 0; 0 1], N[0 2; 1 0]])
    B = MatrixZonotope(N[2 -1; 1 3], [N[0 1; 1 0], N[1 0; 0 -1], N[2 2; -1 0]])

    R = overapproximate(A * B, MatrixZonotope)
    @test center(R) == center(A) * center(B)
    @test ngens(R) == 11   # 3 + 2 + 6
    # spot checks
    @test any(==(center(A) * generators(B)[1]), generators(R))
    @test any(==(generators(A)[2] * center(B)), generators(R))
    @test any(==(generators(A)[1] * generators(B)[3]), generators(R))
end

for N in @tN([Float32, Float64])
    P = SparsePolynomialZonotope(N[1, -1], N[1 1; 0 -1], hcat(N[0, 1]), [2 1; 0 1; 1 0], [1, 3, 5])
    M = N[1 1; -1 1]
    # overapproximate matrix zonotope with interval matrix
    @static if isdefined(@__MODULE__, :IntervalMatrices)
        MZ = MatrixZonotope(N[-1 -4; 4 -1], [N[0.1 0.1; 0.1 0.1]])
        IM = overapproximate(MZ, IntervalMatrix)
        @test IM == IntervalMatrix(N[-1.1 -4.1; 3.9 -1.1], N[-0.9 -3.9; 4.1 -0.9])

        # overapproximate MatrixZonotopeExp
        # test degenerate case
        c = N[1 -2; 2 -1]
        A = MatrixZonotope(c, [zeros(N, 2, 2)])
        expA = MatrixZonotopeExp(A)
        res = overapproximate(expA, MatrixZonotope, 20; tol=1e-5) #for large k it converges to exp(c)
        @test isapprox(center(res), exp(c))
        @test ngens(res) == 0

        # degenerate case + product 
        B = MatrixZonotope(N[1 0; 0 1], Matrix{N}[])
        expAB = MatrixZonotopeExp(A * B)
        res2 = overapproximate(expA, MatrixZonotope, 20)
        @test isapprox(center(res), center(res2))

        # monotonicity test
        @static if isdefined(@__MODULE__, :ExponentialUtilities) || isdefined(@__MODULE__, :Expokit)
            mzexp = MatrixZonotopeExp(MZ)
            em = ExponentialMap(mzexp, P)
            MPex = overapproximate(em, SparsePolynomialZonotope, 2)
            MZex = overapproximate(MPex, Zonotope)

            Z = overapproximate(convert(SimpleSparsePolynomialZonotope, P), Zonotope)
            em_z = ExponentialMap(mzexp, Z)
            Zex = overapproximate(em_z, Zonotope, 5)

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