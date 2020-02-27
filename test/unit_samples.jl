using Distributions: Uniform

P1 = BallInf([0.6, 0.1, -1.3, -0.4, 0.2], 0.6)
A = [1.0  0.0;
     0.0  1.0;
    -1.0  0.0;
     0.0 -1.0]
b = [1.0, 2.0, 3.0, 4.0]
P2 = HPolyhedron(A, b)
P3 = HPolyhedron([0.0 0.1], [3.0])
P4 = Ball2([0.2, -0.3, -1.1, 0.6, -0.7], 0.4)

## Test rand samples are contained in the set
p1 = LazySets.sample(P1)
@test p1 ∈ P1
p1_samples = LazySets.sample(P1, 100)
@test sum(p1_samples .∈ Ref(P1)) == length(p1_samples)

@test LazySets.RejectionSampler(P2).box_approx == [Uniform(-3.0,1.0), Uniform(-4.0,2.0)]
