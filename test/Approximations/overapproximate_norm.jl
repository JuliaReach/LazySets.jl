using LazySets, Test
using LazySets.Approximations: overapproximate_norm

for N in [Float64, Float32]
    for _ in 1:5
        Z = rand(Zonotope; N=N, dim=8, num_generators=5)
        @test overapproximate_norm(Z, 1) ≥ norm(Z, 1)
    end

    Z = Zonotope(N[1, -1], zeros(N, 2, 0))
    @test overapproximate_norm(Z, 1) == 2
end
