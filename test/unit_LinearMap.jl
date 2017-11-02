# π/2 trigonometric rotation
b = BallInf([1., 2.], 1.)
M = [0. -1. ; 1. 0.]
# Test Construction
lm1 = LinearMap(M, b)
@test lm1.M == M
@test lm1.sf == b
# Test Dimension
@test dim(lm1) == 2
# Test Support Vector
d = [1., 1.]
@test σ(d, lm1) == [-1., 2.]
d = [-1., 1.]
@test σ(d, lm1) == [-3., 2.]
d = [-1., -1.]
@test σ(d, lm1) == [-3., 0.]
d = [1., -1.]
@test σ(d, lm1) == [-1., 0.]

# 2D -> 1D Projection
b = BallInf([1., 2.], 1.)
M = [1. 0.]
lm = M*b
# Test Dimension
@test dim(lm) == 1
# Test Support Vector
d = [1.]
@test σ(d, lm) == [2.]
d = [-1.]
@test σ(d, lm) == [0.]

# scalar multiplication
b = BallInf([0., 0.], 1.)
lm = 2.0*b
# Test Dimension
@test dim(lm) == 2
# Test Support Vector
d = [1., 1.]
@test σ(d, lm) == [2., 2.]
d = [-1., 1.]
@test σ(d, lm) == [-2., 2.]
d = [-1., -1.]
@test σ(d, lm) == [-2., -2.]
d = [1., -1.]
@test σ(d, lm) == [2., -2.]

# Nested construction
lm1_copy = LinearMap(eye(2), lm1)
@test lm1_copy.M == lm1.M
@test lm1_copy.sf == lm1.sf
