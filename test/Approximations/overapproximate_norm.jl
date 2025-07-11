using LazySets, Test

for N in [Float64, Float32]
    TOL = 1e-3 #tolerance due to floating-point rounding errors
    for _ in 1:5
        Z = rand(Zonotope; N=N, dim=8, num_generators=5)
        @test overapproximate_norm(Z, 1) ≥ (norm(Z, 1) -TOL)
    end

    Z = Zonotope(N[1, -1], zeros(N, 2, 0))
    @test overapproximate_norm(Z, 1) == 2

    for _ in 1:5
        MZ = rand(MatrixZonotope; N=N, dim=(3,3), num_generators=3)
        @test overapproximate_norm(MZ, 1) ≥ (norm(MZ, 1) -TOL)
    end
end
