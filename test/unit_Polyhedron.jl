# 2D Polyhedron
p = Polyhedron(2)
c1 = LinearConstr([2., 2.], 12.)
c2 = LinearConstr([-3., 3.], 6.)
c3 = LinearConstr([-1., -1.], 0.)
c4 = LinearConstr([2., -4.], 0.)
addconstraint!(p, c3)
addconstraint!(p, c1)
addconstraint!(p, c4)
addconstraint!(p, c2)
# Test Support Vector
d = [1., 0.]
@test σ(d, p) == [4., 2.]
d = [0., 1.]
@test σ(d, p) == [2., 4.]
d = [-1., 0.]
@test σ(d, p) == [-1., 1.]
d = [0., -1.]
@test σ(d, p) == [0., 0.]
