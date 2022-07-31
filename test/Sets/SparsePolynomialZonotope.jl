for N in [Float64, Float32, Rational{Int}]

    # example from Niklas thesis (page 32)
    c = N[4, 4]
    G = N[2 1 2;0 2 2]
    GI = N[1; 0]
    E = [1 0 3;0 1 1]
    PZ = SparsePolynomialZonotope(c, G, GI, E)

    @test center(PZ) == c
    @test dependent_genmat(PZ) == G
    @test independent_genmat(PZ) == GI
    @test expmat(PZ) == E
    @test indexvector(PZ) == [1, 2]

    @test dim(PZ) == 2
    @test ndependentgens(PZ) == 3
    @test nindependentgens(PZ) == 1
    @test ngens(PZ) == 4
    @test nparams(PZ) == 2
    @test order(PZ) == 2//1

    M = N[-0.5 0.2;-0.1 0.6]
    PZ = SparsePolynomialZonotope(zeros(N, 2), N[2 0 1;1 2 1], zeros(N, 2, 0), [1 0 1;0 1 3])
    LMPZ = linear_map(M, PZ)
    @test center(PZ) == zeros(N, 2)
    @test dependent_genmat(LMPZ) ≈ N[-0.8 0.4 -0.3;0.4 1.2 0.5] atol=1e-7
    @test isempty(independent_genmat(LMPZ))
    @test expmat(LMPZ) == [1 0 1;0 1 3]
    @test indexvector(LMPZ) == indexvector(PZ)
end
