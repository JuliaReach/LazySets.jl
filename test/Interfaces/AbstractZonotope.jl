using LazySets, Test
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

for N in @tN([Float64, Float32, Rational{Int}])
    # reduce_order
    c = N[2, 1]
    G = N[-1//2 3//2 1//2 1 0 1; 1//2 3//2 1 1//2 1 0]
    for method in [LazySets.ASB10(), LazySets.COMB03(), LazySets.GIR05(), LazySets.SRMB16()]
        Z = Zonotope(c, N[-1//2 3//2 1//2 1; 1//2 3//2 1 1//2])
        Z2 = reduce_order(Z, 1, method)
        @test ngens(Z2) == 2
        @test order(Z2) == 1
        Z2 = reduce_order(Z, 2, method)
        @test ngens(Z2) == 4
        @test order(Z2) == 2
        Z = Zonotope(c, G)
        Z2 = reduce_order(Z, 2, method)
        @test ngens(Z2) == 4
        @test order(Z2) == 2
        Z2 = Zonotope(N[1, 2], Matrix{N}(undef, 2, 0))
        @test ngens(Z2) == 0
        @test genmat(Z2) == Matrix{N}(undef, 2, 0)
        @test collect(generators(Z2)) == Vector{N}()
    end
    @static if isdefined(@__MODULE__, :StaticArrays)
        using StaticArrays: SVector, SMatrix

        # TODO LazySets.ASB10() does not preserve SMatrix
        # TODO LazySets.SRMB16() should convert to MMatrix
        for method in [LazySets.COMB03(), LazySets.GIR05()]
            Z = Zonotope(SVector{2}(c), SMatrix{2,6}(G))
            @test reduce_order(Z, 2, method) isa Zonotope{N,SVector{2,N},SMatrix{2,4,N,8}}
        end
    end

    # remove_redundant_generators
    M = N[1 0 -2]
    M2 = remove_redundant_generators(M)
    @test M2 isa Matrix{N} && M2 == hcat(N[3])
    M = N[1 0 2 -3 4; 2 0 4 -6 1]
    M2 = remove_redundant_generators(M)
    @test M2 isa Matrix{N} && M2 == N[6 4; 12 1]
end

for N in @tN([Float64, Float32])
    # remove_redundant_generators: optimized
    G = N[1 1 0; 0 0 1] #after sorting and normalization [0 1 1; 1 0 0]
    @test remove_redundant_generators(G) == N[0 2; 1 0]

    H = N[0 0 1; 1 1 0] # equal after sorting
    @test remove_redundant_generators(H) == N[0 1; 2 0]
end
