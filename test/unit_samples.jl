# using Pkg
# Pkg.activate("LazySets")
# using Test

using LazySets
using Distributions
using Polyhedra

n = 5
P1 = rand(BallInf, dim=n)
A = [1.0  0.0;
     0.0  1.0;
    -1.0  0.0;
     0.0 -1.0]
b = [1.0, 2.0, 3.0, 4.0]
P2 = HPolyhedron(A, b)
P3 = HPolyhedron([0.0 0.1], [3.0])
P4 = rand(Ball2, dim=n)

## Test rand samples are contained in the set
p1 = rand(P1)
@test p1 ∈ P1
p1_samples = rand(P1, 100)
@test sum(p1_samples .∈ Ref(P1)) == length(p1_samples)

# test _canonical_length
@test LazySets._canonical_length(P1+-P1.center) ≈ radius(P1)*[-ones(n) ones(n)]
@test LazySets._canonical_length(P2) == [-b[3:4] b[1:2]]
#does not work with box_approximation in _canonical_length at the moment
@test LazySets._canonical_length(P3) == [-Inf Inf; -Inf 30.0]
@test LazySets._canonical_length(P4+-P4.center) ≈ radius(P4)*[-ones(n) ones(n)]
