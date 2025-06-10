using Distributions: Uniform, Normal, MultivariateNormal
using LazySets: DefaultUniform

for N in [Float64]
    P1 = BallInf([0.6, 0.1, -1.3, -0.4, 0.2], 0.6)
    A = [1.0 0.0;
         0.0 1.0;
         -1.0 0.0;
         0.0 -1.0]
    b = [1.0, 2.0, 3.0, 4.0]
    P2 = HPolyhedron(A, b)
    P3 = VPolygon([N[0, -1], N[2, 1]])

    # Test rand samples are contained in the set
    p1 = sample(P1)
    @test p1 ∈ P1
    p1_samples = sample(P1, 100)
    @test sum(p1_samples .∈ Ref(P1)) == length(p1_samples)

    # default distribution
    @test LazySets.RejectionSampler(P2).distribution ==
          [DefaultUniform(-3.0, 1.0), DefaultUniform(-4.0, 2.0)]

    # random-walk sampler
    for P in [P1, P2, P3]
        for sampler in [LazySets.RandomWalkSampler(true), LazySets.RandomWalkSampler(false)]
            p_samples = sample(P, 100; sampler=sampler)
            @test all(v ∈ P for v in p_samples)
        end
    end

    # random-walk sampler
    for P in [P1, P2, P3]
        p_samples = sample(P, 100; sampler=LazySets.CombinedSampler())
        @test all(v ∈ P for v in p_samples)
    end

    # specifying a distribution from Distributions.jl
    @test LazySets.RejectionSampler(P2, Uniform).distribution ==
          [Uniform(-3.0, 1.0), Uniform(-4.0, 2.0)]
    @test LazySets.RejectionSampler(P2, Normal).distribution ==
          [Normal(-3.0, 1.0), Normal(-4.0, 2.0)]

    # test univariate distributions
    X = Interval(1, 2.0)
    v = sample(X, 5)
    @test all(vi ∈ X for vi in v)
    d = Normal(1.5, 0.5)

    v = sample(X, 5; sampler=LazySets.RejectionSampler(d))
    @test all(vi ∈ X for vi in v)

    v = sample(X, 5; sampler=LazySets.RejectionSampler([d]))
    @test all(vi ∈ X for vi in v)

    v = sample(X, 5; sampler=LazySets.RejectionSampler([DefaultUniform(1, 2)]))
    @test all(vi ∈ X for vi in v)

    v = sample(X, 5; sampler=LazySets.RejectionSampler(DefaultUniform(1, 2)))
    @test all(vi ∈ X for vi in v)

    # test multivariate distribution
    H = rand(Hyperrectangle; dim=2)
    sampler = LazySets.RejectionSampler(MultivariateNormal(N[1 0; 0 1.0]); tight=true)
    v = sample(H, 5; sampler=sampler)

    # including vertices
    for k in 0:4
        p1 = sample(P1, 10; include_vertices=k)
        @test length(p1) == 10 + k
    end
    @test length(sample(P1, 10; include_vertices=false)) == 10
    @test length(sample(P1, 10; include_vertices=true)) == 42

    # face sampling
    H = Hyperrectangle(N[1, 2], N[3, 4])
    v0 = sample(H, 10; sampler=LazySets.FaceSampler(0))
    v1 = sample(H, 10; sampler=LazySets.FaceSampler(1))
    v2 = sample(H, 10; sampler=LazySets.FaceSampler(2))
    @test all(all(vi ∈ H for vi in v) for v in [v0, v1, v2])
    vs = vertices_list(H)
    @test all(vi ∈ vs for vi in v0)  # 0-faces are the vertices
    @test all(any(vi[i:i] ∈ project(H, [i]) for i in 1:2) for vi in v1)  # 1-faces are at the border

    # SparsePolynomialZonotope
    c  = N[1.0, 1.0]
    G  = 0.1 * N[2.0 0.0 1.0;
                 1.0 2.0 1.0]
    GI = 0.1 * reshape(N[1.0, 0.5], 2, 1)
    E  = [1 0 1;
           0 1 3]
    idx = [1, 2]
    P4 = SparsePolynomialZonotope(c, G, GI, E, idx)

    # test rand samples
    sampler = LazySets.PolynomialZonotopeSampler()
    pts = sample(P4, 100; sampler = sampler)
    @test length(pts) == 100
    @test_broken all(p ∈ P4 for p in pts)
    P4z = overapproximate(P4, Zonotope)  # temporary sufficient test
    @test all(p ∈ P4z for p in pts)
end
