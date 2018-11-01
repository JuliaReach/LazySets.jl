for _dummy_ in 1:1 # avoid global variable warnings
    # reseeding with random seed
    rng = LazySets.GLOBAL_RNG
    seed = rand(1:10000)
    LazySets.reseed(rng, seed)
    n1 = rand(Int)
    LazySets.reseed(rng, seed)
    n2 = rand(Int)
    @test n1 == n2
end
