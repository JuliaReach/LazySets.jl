using LazySets, Test
using LazySets.Approximations: overapproximate_norm

let
    for _ in 1:5
        Z = rand(Zonotope; dim=8, num_generators=5)
        @test overapproximate_norm(Z, 1) ≥ norm(Z, 1)
    end

    Z = Zonotope([1.0, -1.0], zeros(Float64, 2, 0))
    @test overapproximate_norm(Z, 1) == 2
end
