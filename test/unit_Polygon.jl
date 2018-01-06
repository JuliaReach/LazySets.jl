# Empty polygon
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
@test p.constraints_list[1] == c1
@test p.constraints_list[2] == c2
@test p.constraints_list[3] == c3
@test p.constraints_list[4] == c4

# Test Support Vector
d = [1., 0.]
@test σ(d, p) == [4., 2.]
d = [0., 1.]
@test σ(d, p) == [2., 4.]
d = [-1., 0.]
@test σ(d, p) == [-1., 1.]
d = [0., -1.]
@test σ(d, p) == [0., 0.]

# Test containment
@test ∈([0., 0.], p)
@test ∈([4., 2.], p)
@test ∈([2., 4.], p)
@test ∈([-1., 1.], p)
@test ∈([2., 3.], p)
@test ∈([1., 1.], p)
@test ∈([3., 2.], p)
@test ∈([5./4., 7./4.], p)
@test !∈([4., 1.], p)
@test !∈([5., 2.], p)
@test !∈([3., 4.], p)
@test !∈([-1., 5.], p)

# an_element function
@test an_element(p) ∈ p

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
@test ∈([0., 0.], po)
@test ∈([4., 2.], po)
@test ∈([2., 4.], po)
@test ∈([-1., 1.], po)
@test ∈([2., 3.], po)
@test ∈([1., 1.], po)
@test ∈([3., 2.], po)
@test ∈([5./4., 7./4.], po)
@test !∈([4., 1.], po)
@test !∈([5., 2.], po)
@test !∈([3., 4.], po)
@test !∈([-1., 5.], po)

# an_element function
@test an_element(po) ∈ po

# Test VRepresentation
vp = tovrep(p)
@test [2., 4.] ∈ vp.vertices_list
@test [-1., 1.] ∈ vp.vertices_list
@test [0., 0.] ∈ vp.vertices_list
@test [4., 2.] ∈ vp.vertices_list

# test convex hull of a set of points using the default algorithm
points = [[0.9,0.2], [0.4,0.6], [0.2,0.1], [0.1,0.3], [0.3,0.28]]
vp = VPolygon(points) # by default, a convex hull is run
@test vertices_list(vp) == [ [0.1,0.3],[0.2,0.1], [0.9,0.2],[0.4,0.6] ]

vp = VPolygon(points, apply_convex_hull=false) # we can turn it off
@test vertices_list(vp) == [[0.9,0.2], [0.4,0.6], [0.2,0.1], [0.1,0.3], [0.3,0.28]]

# test for pre-sorted points
points = [[0.1, 0.3], [0.2, 0.1], [0.4, 0.3], [0.4, 0.6], [0.9, 0.2]]
vp = VPolygon(points, algorithm="monotone_chain_sorted")
@test vertices_list(vp) == [[0.1, 0.3], [0.2, 0.1], [0.9, 0.2], [0.4, 0.6]]

# test support vector of a VPolygon
p = HPolygon()
for ci in [c1, c2, c3, c4]
    addconstraint!(p, ci)
end
p = tovrep(p)

# Test Support Vector
d = [1., 0.]
@test σ(d, p) == [4., 2.]
d = [0., 1.]
@test σ(d, p) == [2., 4.]
d = [-1., 0.]
@test σ(d, p) == [-1., 1.]
d = [0., -1.]
@test σ(d, p) == [0., 0.]

# test that #83 is fixed
v = VPolygon([[2.0, 3.0]])
@test [2.0, 3.0] ∈ v
@test !([3.0, 2.0] ∈ v)
v = VPolygon([[2.0, 3.0], [-1.0, -3.4]])
@test [-1.0, -3.4] ∈ v

# an_element function
v = VPolygon([[2., 3.]])
@test an_element(v) ∈ v

# subset
p1 = VPolygon([[0., 0.],[2., 0.]])
p2 = VPolygon([[1., 0.]])
@test ⊆(p2, p1) && ⊆(p2, p1, true)[1]
@test ⊆(HPolygon(), p1)
@test ⊆(p1, BallInf([1., 0.], 1.))
