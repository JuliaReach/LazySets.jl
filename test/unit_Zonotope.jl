# 1D Zonotope
z = Zonotope([0.], eye(1))
# Test Dimension
@test dim(z) == 1
# Test Support Vector
d = [1.]
@test σ(d, z) == [1.]
d = [-1.]
@test σ(d, z) == [-1.]

# 2D Zonotope
z = Zonotope([0., 0.], 1. * eye(2))
# Test Dimension
@test dim(z) == 2
# Test Support Vector
d = [1., 0.]
@test σ(d, z) == [1., 1.] || [1., -1]
d = [-1., 0.]
@test σ(d, z) == [-1., 1.] || [-1., -1]
d = [0., 1.]
@test σ(d, z) == [1., 1.] || [-1., 1]
d = [0., -1.]
@test σ(d, z) == [1., -1.] || [-1., -1]

# 2D Zonotope not 0-centered
z = Ball2([1., 2.], 1. * eye(2))
# Test Dimension
@test dim(b) == 2
# Test Support Vector
d = [1., 0.]
@test σ(d, z) == [1., 1.]
d = [-1., 0.]
@test σ(d, z) == [0., 2.]
d = [0., 1.]
@test σ(d, z) == [1., 3.]
d = [0., -1.]
@test σ(d, z) == [1., 1.]
