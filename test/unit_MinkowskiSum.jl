# Sum of 2D centered balls in norm 2 and infinity
b1 = BallInf([0., 0.], 2.)
b2 = Ball2([0., 0.], 1.)
# Test Construction
X = MinkowskiSum(b1, b2)
@test X.X == b1
@test X.Y == b2
# Test Dimension
@test dim(X) == 2
# Test Support Vector
d = [1., 0.]
v = σ(d, X)
@test v[1] == 3.
d = [-1., 0.]
v = σ(d, X)
@test v[1] == -3.
d = [0., 1.]
v = σ(d, X)
@test v[2] == 3.
d = [0., -1.]
v = σ(d, X)
@test v[2] == -3.

# Sum of not-centered 2D balls in norm 2 and infinity
b1 = BallInf([-1., 3.], 2.)
b2 = Ball2([1., 2.], 1.)
s = b1 + b2
# Test Support Vector
d = [1., 0.]
v = σ(d, s)
@test v[1] == 3.
d = [-1., 0.]
v = σ(d, s)
@test v[1] == -3.
d = [0., 1.]
v = σ(d, s)
@test v[2] == 8.
d = [0., -1.]
v = σ(d, s)
@test v[2] == 2.

# Sum of array of LazySet
# 2-elements
ms = MinkowskiSum(Singleton([1.]), Singleton([2.]))
@test ρ([1.], ms) == 3.
@test ρ([-1.], ms) == -3.
# 3-elements
ms = MinkowskiSum(Singleton([1.]), MinkowskiSum(Singleton([2.]), Singleton([3.])))
@test ρ([1.], ms) == 6.
@test ρ([-1.], ms) == -6.
