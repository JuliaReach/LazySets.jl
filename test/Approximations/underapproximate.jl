for N in [Float64, Rational{Int}, Float32]
    X = Hyperrectangle(zeros(N, 2), ones(N, 2))
    U = underapproximate(X, BoxDirections{N}(2))
    @test U ⊆ X

    X = VPolygon([N[4, 2], N[2, 0], N[0, 2], N[2, 4]])
    U = underapproximate(X, Hyperrectangle)
    @test U ≈ Hyperrectangle(N[2, 2], ones(N, 2))
end

for N in [Float64, Float32]
    X = VPolygon([N[1, 0], N[1, 2], N[-1, 2], N[-1, 1//3], N[-2//3, 0]])
    U = underapproximate(X, Ball2)
    @test U ≈ Ball2(N[0, 1], N(1))
end
