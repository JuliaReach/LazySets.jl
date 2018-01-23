for N in [Float64, Rational{Int}, Float32]
    B = BallInf(ones(N, 2), N(3.))
    H = Hyperrectangle(ones(N, 2), ones(N, 2))
    I = Intersection(B, H)

    # dim
    @test dim(I) == 2

    # support vector (error)
    @test_throws ErrorException σ(ones(N, 2), I)

    # membership
    @test ∈(ones(N, 2), I) && !∈(N[5., 5.], I)

    # emptiness of intersection
    @test !isempty(I)
end
