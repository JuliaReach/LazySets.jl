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
z = Zonotope([1., 2.], 1. * eye(2))
# Test Dimension
@test dim(b) == 2
# Test Support Vector
d = [1., 0.]
@test σ(d, z) == [2., 3]
d = [-1., 0.]
@test σ(d, z) == [0., 3]
d = [0., 1.]
@test σ(d, z) == [2., 3]
d = [0., -1.]
@test σ(d, z) == [2., 1.]

# an_element function
@test an_element(z) ∈ z

# concrete operations
Z1 = Zonotope([1, 1.], [1 1; -1 1.])
Z2 = Zonotope([-1, 1.], eye(2))
A = [0.5 1; 1 0.5]

# concrete Minkowski sum
Z3 = minkowski_sum(Z1, Z2)
@test Z3.center == [0, 2.]
@test Z3.generators == [1 1 1 0; -1 1 0 1.]

# concrete linear map and scale
Z4 = linear_map(A, Z3)
@test Z4.center == [2, 1.]
@test Z4.generators == [-0.5 1.5 0.5 1.0; 0.5 1.5 1.0 0.5]
Z5 = scale(0.5, Z3)
@test Z5.center == [0, 1.]
@test Z5.generators == [0.5 0.5 0.5 0.0; -0.5 0.5 0.0 0.5]
