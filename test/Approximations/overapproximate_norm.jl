using LazySets, Test
using LazySets.ReachabilityBase.Comparison: _geq

for N in [Float64, Float32]
    for _ in 1:5
        Z = rand(Zonotope; N=N, dim=8, num_generators=5)
        @test _geq(overapproximate_norm(Z, 1), norm(Z, 1), atol=1e-3)
    end

    Z = Zonotope(N[1, -1], zeros(N, 2, 0))
    @test overapproximate_norm(Z, 1) == 2

    for _ in 1:5
        MZ = rand(MatrixZonotope; N=N, dim=(3, 3), num_generators=3)
        @test _geq(overapproximate_norm(Z, 1), norm(Z, 1), atol=1e-3)
    end
end
