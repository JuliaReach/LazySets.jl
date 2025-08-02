using LazySets, Test
if !isdefined(@__MODULE__, Symbol("@tN"))
    macro tN(v)
        return v
    end
end

for N in @tN([Float64, Float32, Rational{Int}])
    # π/2 trigonometric rotation
    b = BallInf(N[1, 2], N(1))
    M = N[0 -1; 1 0]

    # construction
    lm1 = InverseLinearMap(M, b)
    @test lm1.M == M
    @test lm1.X == b

    # dimension
    @test dim(lm1) == 2

    # Test Support Vector and Support Function
    d = N[1, 1]
    @test σ(d, lm1) == N[3, 0]
    @test ρ(d, lm1) == N(3)
    d = N[-1, 1]
    @test σ(d, lm1) == N[1, 0]
    @test ρ(d, lm1) == N(-1)
    d = N[-1, -1]
    @test σ(d, lm1) == N[1, -2]
    @test ρ(d, lm1) == N(1)
    d = N[1, -1]
    @test σ(d, lm1) == N[3, -2]
    @test ρ(d, lm1) == N(5)

    # scalar multiplication and concretize
    b = BallInf(N[0, 0, 0], N(2))
    ilm1 = InverseLinearMap(N(5), b)
    ilm2 = InverseLinearMap(Matrix(ilm1.M), b)
    @test ilm1 == ilm2
    @test overapproximate(concretize(ilm2.M * ilm2), Hyperrectangle) ==
          overapproximate(b, Hyperrectangle)

    # emptiness check
    @test !isempty(ilm1)
end
