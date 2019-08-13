for N in [Float64, Rational{Int}, Float32]
    X = Hyperrectangle(zeros(N, 2), ones(N, 2))
    U = underapproximate(X, BoxDirections{N}(2))
    @test U ⊆ X
end
