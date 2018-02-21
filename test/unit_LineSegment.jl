for N in [Float64, Rational{Int}, Float32]
    # construction
    p, q = N[1., 1.], N[2., 2.]
    l = LineSegment(p, q)

    # dimension
    @test dim(l) == 2

    # support vector
    @test σ(N[1., 1.], l) == q
    @test σ(N[0., 1.], l) == q
    @test σ(N[-1., 1.], l) == q
    @test σ(N[-1., 0.], l) == p
    @test σ(N[-1., -1.], l) == p
    @test σ(N[0., -1.], l) == p
    @test σ(N[1., -1.], l) == q
    @test σ(N[1., 0.], l) == q

    # membership
    @test !∈(N[0., 0.], l)
    @test ∈(N[1., 1.], l)
    @test ∈(N[1.5, 1.5], l)
    @test !∈(N[3., 4.], l)
    @test !∈(N[7., 4.], l)
    @test !∈(N[1.5, 1.6], l)
end
