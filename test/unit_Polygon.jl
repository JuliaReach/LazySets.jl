# Polygon
p = HPolygon()
# Test Dimension
@test dim(p) == 2
# Test that constraints are automatically sorted
c1 = LinearConstraint([2., 2.], 12.)
c2 = LinearConstraint([-3., 3.], 6.)
c3 = LinearConstraint([-1., -1.], 0.)
c4 = LinearConstraint([2., -4.], 0.)
addconstraint!(p, c3)
addconstraint!(p, c1)
addconstraint!(p, c4)
addconstraint!(p, c2)
@test p.constraints[1] == c1
@test p.constraints[2] == c2
@test p.constraints[3] == c3
@test p.constraints[4] == c4
# Test Support Vector
d = [1., 0.]
@test σ(d, p) == [4., 2.]
d = [0., 1.]
@test σ(d, p) == [2., 4.]
d = [-1., 0.]
@test σ(d, p) == [-1., 1.]
d = [0., -1.]
@test σ(d, p) == [0., 0.]
# Test VRepresentation
vp = tovrep(p)
@test vp.vl[1] == [2., 4.]
@test vp.vl[2] == [-1., 1.]
@test vp.vl[3] == [0., 0.]
@test vp.vl[4] == [4., 2.]

# Test containment
@test is_contained([0., 0.], p)
@test is_contained([4., 2.], p)
@test is_contained([2., 4.], p)
@test is_contained([-1., 1.], p)
@test is_contained([2., 3.], p)
@test is_contained([1., 1.], p)
@test is_contained([3., 2.], p)
@test is_contained([5./4., 7./4.], p)
@test !is_contained([4., 1.], p)
@test !is_contained([5., 2.], p)
@test !is_contained([3., 4.], p)
@test !is_contained([-1., 5.], p)

# Optimized Polygon
po = HPolygonOpt(p)
# Test Dimension
@test dim(po) == 2
# Test Support Vector
d = [1., 0.]
@test σ(d, po) == [4., 2.]
d = [0., 1.]
@test σ(d, po) == [2., 4.]
d = [-1., 0.]
@test σ(d, po) == [-1., 1.]
d = [0., -1.]
@test σ(d, po) == [0., 0.]
# Test containment

@test is_contained([0., 0.], po)
@test is_contained([4., 2.], po)
@test is_contained([2., 4.], po)
@test is_contained([-1., 1.], po)
@test is_contained([2., 3.], po)
@test is_contained([1., 1.], po)
@test is_contained([3., 2.], po)
@test is_contained([5./4., 7./4.], po)
@test !is_contained([4., 1.], po)
@test !is_contained([5., 2.], po)
@test !is_contained([3., 4.], po)
@test !is_contained([-1., 5.], po)

# test conveh hull algorithm using andrew_monotone_chain
points = [[0.9,0.2], [0.4,0.6], [0.2,0.1], [0.1,0.3], [0.3,0.28]]
@test andrew_monotone_chain(points) == [ [0.1,0.3],[0.2,0.1], [0.9,0.2],[0.4,0.6] ]
