using LazySets, Test

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

    # transpose
    MZt = transpose(MZ)
    @test center(MZt) == c
    @test generators(MZt) == [transpose(Ai) for Ai in generators(MZ)]

    #norm
    @test norm(MZ, Inf) == 6
    @test norm(MZ, 1) == 6
    
    #linear_map 
    M = N[2 -1; 1 0]
    lm_l = linear_map(M, MZ)
    @test center(lm_l) == [2 -1; 1 0]
    @test generators(lm_l) == [[2 -4; 1 -1], [5 -3; 2 -1]]
    lm_r = linear_map(MZ, M)
    @test center(lm_r) == [2 -1; 1 0]
    @test generators(lm_r) == [[1 -1; 2 0], [3 -2; -1 1]]
end

for N in (Float64, Float32, Rational{Int})
    # A: (2x2), B: (2x3)
    A_c = N[1 0; 0 1]
    A_gens = [N[1 -1; 0 2], N[2 -1; -1 1]]
    A_mz = MatrixZonotope(A_c, A_gens)

    B_c = N[1 0 0; 0 1 1]
    B_gens = [N[1 0 -1; -1 2 0], N[1 1 -1; 0 2 -1]]
    B_mz = MatrixZonotope(B_c, B_gens)

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
    MZP_chain = A_mz * B_mz * MatrixZonotope(N[1 0; 0 1; 1 1], [N[1 1; 0 0; 0 -2], N[0 1; 1 0; 1 -1]])
    @test nfactors(MZP_chain) == 3
    @test MZP_chain isa MatrixZonotopeProduct{N}

    # remove_redundant_factors
    const_A = MatrixZonotope(N[1 0; 0 1], Vector{Matrix{N}}())
    mixed = MatrixZonotopeProduct(const_A, A_mz, B_mz)
    simplified = remove_redundant_factors(mixed)
    @test nfactors(simplified) == 2
    @test size(simplified) == (2, 3)

    # incompatible dimensions
    C_c = N[1 2; 3 4; 5 6]
    C_gens = [N[2 0; 1 -1; -3 0]]
    C_mz = MatrixZonotope(C_c, C_gens)  # 3x2
    @test_throws AssertionError MatrixZonotopeProduct(A_mz, C_mz)
end

for N in @tN([Float64, Float32])
    MZ = rand(MatrixZonotope; N=N)
    @test MZ isa MatrixZonotope{N}
end
