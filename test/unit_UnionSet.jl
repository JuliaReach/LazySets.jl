# UnionSet of Polygon
# Union construction
p1 = HPolygon()
c = LinearConstr([1., 0.], 1.)
addconstraint!(p1, c)
c = LinearConstr([-1., 0.], 1.)
addconstraint!(p1, c)
c = LinearConstr([0., 1.], 1.)
addconstraint!(p1, c)
c = LinearConstr([0., -1.], 1.)
addconstraint!(p1, c)
p2 = HPolygon()
v = [4., 5.]
c = LinearConstr([1., 0.], 1. + dot(v, [1., 0.]))
addconstraint!(p2, c)
c = LinearConstr([-1., 0.], 1. + dot(v, [-1., 0.]))
addconstraint!(p2, c)
c = LinearConstr([0., 1.], 1. + dot(v, [0., 1.]))
addconstraint!(p2, c)
c = LinearConstr([0., -1.], 1. + dot(v, [0., -1.]))
addconstraint!(p2, c)
u = UnionSet([p1, p2])
# In test
@test in([0., 0.], u) == true
@test in([1., 0.], u) == true
@test in([-1., 0.], u) == true
@test in([0., 1.], u) == true
@test in([0., -1.], u) == true
@test in([4., 5.], u) == true
@test in([5., 5.], u) == true
@test in([3., 5.], u) == true
@test in([4., 6.], u) == true
@test in([4., 4.], u) == true
@test in([2., 2.], u) == false
@test in([2., 5.], u) == false
@test in([5., 2.], u) == false
@test in([-1., 2.], u) == false
@test in([-5., -5.], u) == false
