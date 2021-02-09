using Distributions: Uniform, Normal, MultivariateNormal
using LazySets: DefaultUniform

for N in [Float64]
    P1 = BallInf([0.6, 0.1, -1.3, -0.4, 0.2], 0.6)
    A = [1.0  0.0;
         0.0  1.0;
        -1.0  0.0;
         0.0 -1.0]
    b = [1.0, 2.0, 3.0, 4.0]
    P2 = HPolyhedron(A, b)
    P3 = VPolygon([N[0, -1], N[2, 1]])

    # Test rand samples are contained in the set
    p1 = LazySets.sample(P1)
    @test p1 ∈ P1
    p1_samples = LazySets.sample(P1, 100)
    @test sum(p1_samples .∈ Ref(P1)) == length(p1_samples)

    # default distribution
    @test LazySets.RejectionSampler(P2).distribution == [DefaultUniform(-3.0,1.0), DefaultUniform(-4.0,2.0)]

    # random-walk sampler
    for P in [P1, P2, P3]
        for sampler in [LazySets.RandomWalkSampler(true), LazySets.RandomWalkSampler(false)]
            p_samples = LazySets.sample(P, 100, sampler=sampler)
            @test all(v ∈ P for v in p_samples)
        end
    end

    # specifying a distribution from Distributions.jl
    @test LazySets.RejectionSampler(P2, Uniform).distribution == [Uniform(-3.0, 1.0), Uniform(-4.0,2.0)]
    @test LazySets.RejectionSampler(P2, Normal).distribution == [Normal(-3.0, 1.0), Normal(-4.0, 2.0)]

    # test univariate distributions
    X = Interval(1, 2.)
    v = sample(X, 5)
    @test all(vi ∈ X for vi in v)
    d = Normal(1.5, 0.5)

    v = sample(X, 5, sampler=LazySets.RejectionSampler(d))
    @test all(vi ∈ X for vi in v)

    v = sample(X, 5, sampler=LazySets.RejectionSampler([d]))
    @test all(vi ∈ X for vi in v)

    v = sample(X, 5, sampler=LazySets.RejectionSampler([DefaultUniform(1, 2)]))
    @test all(vi ∈ X for vi in v)

    v = sample(X, 5, sampler=LazySets.RejectionSampler(DefaultUniform(1, 2)))
    @test all(vi ∈ X for vi in v)

    # test multivariate distribution
    H = rand(Hyperrectangle, dim=2)
    sampler = LazySets.RejectionSampler(MultivariateNormal(N[1 0; 0 1.]), tight=true)
    v = sample(H, 5, sampler=sampler)

    # including vertices
    for k in 0:4
        p1 = sample(P1, 10; include_vertices=k)
        @test length(p1) == 10 + k
    end
    @test length(sample(P1, 10; include_vertices=false)) == 10
    @test length(sample(P1, 10; include_vertices=true)) == 42
end
