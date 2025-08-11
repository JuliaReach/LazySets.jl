using LazySets, Test
using LazySets.MatrixZonotopeModule: vectorize, matrixize
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

for N in @tN([Float64, Float32, Rational{Int}])
    # test constructor
    c = N[1 0; 0 3]
    gens = [N[1 -1; 0 2], N[2 -1; -1 1]]
    MZ = MatrixZonotope(c, gens)
    @test MZ isa MatrixZonotope{N}
    @test center(MZ) == c
    @test generators(MZ) == gens
    @test indexvector(MZ) == [1, 2]
    @test ngens(MZ) == 2

    # test constructor with custom indexvector
    idx = [10, 20]
    MZ2 = MatrixZonotope(c, gens, idx)
    @test indexvector(MZ2) == idx

    @test size(MZ) == size(c)
    @test eltype(MZ) == N

    @test_throws ArgumentError MatrixZonotope(c, gens, [1, 1])
    @test_throws ArgumentError MatrixZonotope(c, gens, [1])

    # test for mismatching matrix sizes
    c_bad_size = N[1 0; 0 1; 0 0]
    @test_throws ArgumentError MatrixZonotope(c_bad_size, gens)

    gens_bad_size = [N[1 -1; 0 2], N[2 -1; -1 1; 0 0]]
    @test_throws ArgumentError MatrixZonotope(c, gens_bad_size)

    # test with empty generators
    gens = Vector{Matrix{N}}()
    MZ3 = MatrixZonotope(c, gens)
    @test MZ3 isa MatrixZonotope{N}
    @test center(MZ3) == c
    @test isempty(generators(MZ3))
    @test isempty(indexvector(MZ3))

    #scale
    MZ_scaled = scale(2, MZ)
    @test MZ_scaled isa MatrixZonotope{N}
    @test center(MZ_scaled) == 2 .* center(MZ)
    @test generators(MZ_scaled) == [2 .* g for g in generators(MZ)]
    # scale!
    MZ_copy = copy(MZ)
    scale!(3, MZ_copy)
    @test center(MZ_copy) == 3 .* center(MZ)
    @test generators(MZ_copy) == [3 .* g for g in generators(MZ)]
    # case: no generators
    MZ4 = scale(2, MZ3)
    @test center(MZ4) == 2 .* c
    @test isempty(generators(MZ4))

    # transpose
    MZt = transpose(MZ)
    @test center(MZt) == c
    @test generators(MZt) == [transpose(Ai) for Ai in generators(MZ)]

    # norm
    @test norm(MZ, Inf) == 7
    @test norm(MZ, 1) == 8

    # linear_map
    M = N[2 -1; 1 0]
    lm_l = linear_map(M, MZ)
    @test center(lm_l) == [2 -3; 1 0]
    @test generators(lm_l) == [[2 -4; 1 -1], [5 -3; 2 -1]]
    lm_r = linear_map(MZ, M)
    @test center(lm_r) == [2 -1; 3 0]
    @test generators(lm_r) == [[1 -1; 2 0], [3 -2; -1 1]]
    @test indexvector(lm_l) == indexvector(lm_r) && indexvector(lm_r) == indexvector(MZ)

    # A: (2x2), B: (2x3)
    A_c = N[1 0; 0 1]
    A_gens = [N[1 -1; 0 2], N[2 -1; -1 1]]
    A_mz = MatrixZonotope(A_c, A_gens)

    # B: (2x3)
    B_c = N[1 0 0; 0 1 1]
    B_gens = [N[1 0 -1; -1 2 0], N[1 1 -1; 0 2 -1]]
    B_mz = MatrixZonotope(B_c, B_gens)

    # C: (3, 2)
    # incompatible dimensions
    C_c = N[1 2; 3 4; 5 6]
    C_gens = [N[2 0; 1 -1; -3 0]]
    C_mz = MatrixZonotope(C_c, C_gens)
    @test_throws AssertionError MatrixZonotopeProduct(A_mz, C_mz)

    # constructor from 2
    MZP = MatrixZonotopeProduct(A_mz, B_mz)
    @test MZP isa MatrixZonotopeProduct{N}
    @test nfactors(MZP) == 2
    @test factors(MZP) == [A_mz, B_mz]

    # size checks
    @test size(MZP) == (2, 3)
    @test size(MZP, 1) == 2
    @test size(MZP, 2) == 3

    # chainable multiplication
    MZP_chain = A_mz * B_mz *
                MatrixZonotope(N[1 0; 0 1; 1 1], [N[1 1; 0 0; 0 -2], N[0 1; 1 0; 1 -1]])
    MZP2 = C_mz * B_mz
    @test nfactors(MZP_chain) == 3
    @test MZP_chain isa MatrixZonotopeProduct{N}
    @test A_mz * B_mz == MZP
    @test A_mz * MZP == MatrixZonotopeProduct(A_mz, A_mz, B_mz)
    @test MZP * C_mz == MatrixZonotopeProduct(A_mz, B_mz, C_mz)
    @test MZP * MZP2 == MatrixZonotopeProduct(A_mz, B_mz, C_mz, B_mz)

    # MatrixZonotopeExp constructor
    expMZ = MatrixZonotopeExp(A_mz)
    @test expMZ.M == A_mz
    @test size(expMZ) == size(A_mz)
    @test_throws AssertionError MatrixZonotopeExp(B_mz)

    # vectorize and matrixize
    c = N[1 0; 0 3]
    gens = [N[1 -1; 0 2]]
    MZ = MatrixZonotope(c, gens)
    Z = vectorize(MZ)
    @test Z == Zonotope(N[1, 0, 0, 3], hcat(N[1, 0, -1, 2]))
    MZ2 = matrixize(Z, (2, 2))
    @test isapprox(center(MZ), center(MZ2))
    @test isapprox(generators(MZ), generators(MZ2))

    # order
    @test order(MZ) == 1//4 && order(MZ) == order(Z)

    # remove redundant generators 
    MZ2 = MatrixZonotope(c, [N[1 4; 0 -2], N[-1 1; 0 -1], N[1 -1; 0 1]])
    MZred = remove_redundant_generators(MZ2)
    @test ngens(MZred) == 2
    # `remove_redundant_generators`` introduces floating point errors for Float32
    # and `isequivalent` is not robust to minor imprecisions 
    if N==Float64 
        @static if isdefined(@__MODULE__, :Polyhedra)
            @test isequivalent(vectorize(MZ2), vectorize(MZred))
        end
    end

    # minkowski sum
    ms =minkowski_sum(MZ, MZ2) 
    @test center(ms) == N[2 0; 0 6]
    @test generators(ms) == [N[1 -1; 0 2], N[1 4; 0 -2], N[-1 1; 0 -1], N[1 -1; 0 1]]
end

for N in @tN([Float64, Float32])
    MZ = rand(MatrixZonotope, N=N, dim=(2,2), num_generators=8)
    @test MZ isa MatrixZonotope{N}

    # reduce order
    MZred = reduce_order(MZ, 1)
    @test order(MZred) ≤ 1
    @test center(MZred) == center(MZ)
end
