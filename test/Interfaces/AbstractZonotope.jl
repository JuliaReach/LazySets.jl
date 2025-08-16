using LazySets, Test
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

for N in @tN([Float64, Float32, Rational{Int}])
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
