# 2D polytope
p = HPolytope()
c1 = LinearConstraint([2., 2.], 12.)
c2 = LinearConstraint([-3., 3.], 6.)
c3 = LinearConstraint([-1., -1.], 0.)
c4 = LinearConstraint([2., -4.], 0.)
addconstraint!(p, c3)
addconstraint!(p, c1)
addconstraint!(p, c4)
addconstraint!(p, c2)

# support vector
d = [1., 0.]
@test σ(d, p) == [4., 2.]
d = [0., 1.]
@test σ(d, p) == [2., 4.]
d = [-1., 0.]
@test σ(d, p) == [-1., 1.]
d = [0., -1.]
@test σ(d, p) == [0., 0.]

# membership
@test ∈([5./4., 7./4.], p)
@test !∈([4., 1.], p)
