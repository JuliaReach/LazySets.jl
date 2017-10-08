# 1D BallInf
b = BallInf([0.], 1.)
# Test Dimension
@test dim(b) == 1
# Test Support Vector
d = [1.]
@test σ(d, b) == [1.]
d = [-1.]
@test σ(d, b) == [-1.]

# 2D BallInf
b = BallInf([0., 0.], 1.)
# Test Dimension
@test dim(b) == 2
# Test Support Vector
d = [1., 1.]
@test σ(d, b) == [1., 1.]
d = [-1., 1.]
@test σ(d, b) == [-1., 1.]
d = [-1., -1.]
@test σ(d, b) == [-1., -1.]
d = [1., -1.]
@test σ(d, b) == [1., -1.]

# 2D BallInf not 0-centered
b = BallInf([1., 2.], 1.)
# Test Dimension
@test dim(b) == 2
# Test Support Vector
d = [1., 1.]
@test σ(d, b) == [2., 3.]
d = [-1., 1.]
@test σ(d, b) == [0., 3.]
d = [-1., -1.]
@test σ(d, b) == [0., 1.]
d = [0., -1.]
@test σ(d, b) == [2., 1.]

# 2D BallInf radius =/= 1
b = BallInf([0., 0.], 2.)
# Test Dimension
@test dim(b) == 2
# Test Support Vector
d = [1., 1.]
@test σ(d, b) == [2., 2.]
d = [-1., 1.]
@test σ(d, b) == [-2., 2.]
d = [-1., -1.]
@test σ(d, b) == [-2., -2.]
d = [1., -1.]
@test σ(d, b) == [2., -2.]
