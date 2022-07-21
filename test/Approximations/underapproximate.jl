for N in [Float64, Rational{Int}, Float32]
    X = Hyperrectangle(zeros(N, 2), ones(N, 2))
    U = underapproximate(X, BoxDirections{N}(2))
    @test U ⊆ X

    X = VPolygon([N[4, 2], N[2, 0], N[0, 2], N[2, 4]])
    U = underapproximate(X, Hyperrectangle)
    @test U ≈ Hyperrectangle(N[2, 2], ones(N, 2))
end
