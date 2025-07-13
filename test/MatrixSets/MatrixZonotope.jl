using LazySets, Test

for N in [Float64, Float32, Rational{Int}]
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
    @test generators(MZt)== [transpose(Ai) for Ai in generators(MZ)]

    # norm
    @test norm(MZ, Inf) == 7
    @test norm(MZ, 1) == 8
end

for N in [Float64, Float32]
    MZ = rand(MatrixZonotope; N=N)
    @test MZ isa MatrixZonotope{N}
end
