# Cartesian Product of a centered 1D BallInf and a centered 2D BallInf
# Here a 3D BallInf
b1 = BallInf([0.], 1.)
b2 = BallInf([0., 0.], 1.)
# Test Construction
c1 = CartesianProduct(b1, b2)
@test c1.X == b1
@test c1.Y == b2
# Test Dimension
@test dim(c1) == 3
# Test Support Vector
d = [1., 1., 1.]
@test σ(d, c1) == [1., 1., 1.]
d = [1., 1., -1.]
@test σ(d, c1) == [1., 1., -1.]
d = [1., -1., 1.]
@test σ(d, c1) == [1., -1., 1.]
d = [1., -1., -1.]
@test σ(d, c1) == [1., -1., -1.]
d = [-1., 1., 1.]
@test σ(d, c1) == [-1., 1., 1.]
d = [-1., 1., -1.]
@test σ(d, c1) == [-1., 1., -1.]
d = [-1., -1., 1.]
@test σ(d, c1) == [-1., -1., 1.]
d = [-1., -1., -1.]
@test σ(d, c1) == [-1., -1., -1.]

# Cartesian Product of a not-centered 1D BallInf and a not-centered 2D BallInf
# Here a  Hyperrectangle where c = [1., -3., 4.] and r = [3., 2., 2.]
b1 = BallInf([1.], 3.)
b2 = BallInf([-3., 4.], 2.)
# Test Construction
c = CartesianProduct(b1, b2)
# Test Dimension
@test dim(c) == 3
# Test Support Vector
d = [1., 1., 1.]
@test σ(d, c) == [4., -1., 6.]
d = [1., 1., -1.]
@test σ(d, c) == [4., -1., 2.]
d = [1., -1., 1.]
@test σ(d, c) == [4., -5., 6.]
d = [1., -1., -1.]
@test σ(d, c) == [4., -5., 2.]
d = [-1., 1., 1.]
@test σ(d, c) == [-2., -1., 6.]
d = [-1., 1., -1.]
@test σ(d, c) == [-2., -1., 2.]
d = [-1., -1., 1.]
@test σ(d, c) == [-2., -5., 6.]
d = [-1., -1., -1.]
@test σ(d, c) == [-2., -5., 2.]

# Test Cartesian Product with VoidSet
s = Singleton([1.])
cs1 = VoidSet(1)*s
cs2 = s*VoidSet(1)
@test cs1 isa VoidSet
@test cs2 isa VoidSet

# Test Cartesian Product of an array
# 0-elements
as = LazySet[]
cs = CartesianProduct(as)
@test cs isa VoidSet
# 1-element
as = LazySet[Singleton([1.])]
cs = CartesianProduct(as)
@test cs.element == [1.]
# 3-elements
as = LazySet[Singleton([1.]), Singleton([2.]), Singleton([3.])]
cs = CartesianProduct(as)
@test cs.X.element == [1.]
@test cs.Y.X.element == [2.]
@test cs.Y.Y.element == [3.]

# Test containment with respect to CartesianProduct
p1 = HPolygon()
addconstraint!(p1, LinearConstraint([2., 2.], 12.))
addconstraint!(p1, LinearConstraint([-3., 3.], 6.))
addconstraint!(p1, LinearConstraint([-1., -1.], 0.))
addconstraint!(p1, LinearConstraint([2., -4.], 0.))
p2 = HPolygon()
addconstraint!(p2, LinearConstraint([1., 0.], 1.))
addconstraint!(p2, LinearConstraint([-1., 0.], 1.))
addconstraint!(p2, LinearConstraint([0., 1.], 1.))
addconstraint!(p2, LinearConstraint([0., -1.], 1.))
cp = CartesianProduct(p1, p2)
@test is_contained([0., 0., 0., 0.], cp)
@test is_contained([4., 2., 1., 0.], cp)
@test is_contained([2., 4., -1., 0.], cp)
@test is_contained([-1., 1., .5, .7], cp)
@test is_contained([2., 3., -.8, .9], cp)
@test is_contained([1., 1., -1., 0.], cp)
@test is_contained([3., 2., 0., 1.], cp)
@test is_contained([5./4., 7./4., 1., 1.], cp)
@test !is_contained([4., 1., 0., 0.], cp)
@test !is_contained([5., 2., 0., 0.], cp)
@test !is_contained([3., 4., 0., 0.], cp)
@test !is_contained([-1., 5., 0., 0.], cp)
@test !is_contained([4., 2., 3., 1.], cp)
@test !is_contained([2., 3., 3., 1.], cp)
@test !is_contained([0., 0., 3., 1.], cp)
@test !is_contained([1., 1., 3., 1.], cp)

