using Distributions: Uniform, Normal
using LazySets: DefaultUniform

for N in [Float64]
    P1 = BallInf([0.6, 0.1, -1.3, -0.4, 0.2], 0.6)
    A = [1.0  0.0;
         0.0  1.0;
        -1.0  0.0;
         0.0 -1.0]
    b = [1.0, 2.0, 3.0, 4.0]
    P2 = HPolyhedron(A, b)
    P3 = HPolyhedron([0.0 0.1], [3.0])
    P4 = Ball2([0.2, -0.3, -1.1, 0.6, -0.7], 0.4)

    # Test rand samples are contained in the set
    p1 = LazySets.sample(P1)
    @test p1 ∈ P1
    p1_samples = LazySets.sample(P1, 100)
    @test sum(p1_samples .∈ Ref(P1)) == length(p1_samples)

    # default distribution
    @test LazySets.RejectionSampler(P2).box_approx == [DefaultUniform(-3.0,1.0), DefaultUniform(-4.0,2.0)]

    # polytope sampler
    p1_samples = LazySets.sample(P1, 100, sampler=LazySets.PolytopeSampler)
    @test all(v ∈ P1 for v in p1_samples)

    # specifying a distribution from Distributions.jl
    @test LazySets.RejectionSampler(P2, Uniform).box_approx == [Uniform(-3.0, 1.0), Uniform(-4.0,2.0)]
    @test LazySets.RejectionSampler(P2, Normal).box_approx == [Normal(-3.0, 1.0), Normal(-4.0, 2.0)]

    # including vertices
    for k in 0:4
        p1 = sample(P1, 10; include_vertices=k)
        @test length(p1) == 10 + k
    end
    @test length(sample(P1, 10; include_vertices=false)) == 10
    @test length(sample(P1, 10; include_vertices=true)) == 42
end
