for _dummy_ in 1:1 # avoid global variable warnings
    # reseeding with random seed
    rng = LazySets.GLOBAL_RNG
    seed = rand(1:10000)
    LazySets.reseed(rng, seed)
    n1 = rand(Int)
    LazySets.reseed(rng, seed)
    n2 = rand(Int)
    @test n1 == n2

    # StrictlyIncreasingIndices
    vectors = Vector{AbstractVector{Int}}()
    for v in LazySets.StrictlyIncreasingIndices(5, 4)
        push!(vectors, copy(v))
    end
    @test vectors == [[1, 2, 3, 4], [1, 2, 3, 5], [1, 2, 4, 5], [1, 3, 4, 5], [2, 3, 4, 5]]

    # invertible matrix
    @test LazySets.isinvertible([2 3; 1 2])
    @test !LazySets.isinvertible([2 3; 0 0])
end
