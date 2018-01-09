# 1D Ball2
b = Ball2([0.], 1.)
# Test Dimension
@test dim(b) == 1
# Test Support Vector
d = [1.]
@test σ(d, b) == [1.]
d = [-1.]
@test σ(d, b) == [-1.]

# 2D Ball2
b = Ball2([0., 0.], 1.)
# Test Dimension
@test dim(b) == 2
# Test Support Vector
d = [1., 0.]
@test σ(d, b) == [1., 0.]
d = [-1., 0.]
@test σ(d, b) == [-1., 0.]
d = [0., 1.]
@test σ(d, b) == [0., 1.]
d = [0., -1.]
@test σ(d, b) == [0., -1.]

# 2D Ball2 not 0-centered
b = Ball2([1., 2.], 1.)
# Test Dimension
@test dim(b) == 2
# Test Support Vector
d = [1., 0.]
@test σ(d, b) == [2., 2.]
d = [-1., 0.]
@test σ(d, b) == [0., 2.]
d = [0., 1.]
@test σ(d, b) == [1., 3.]
d = [0., -1.]
@test σ(d, b) == [1., 1.]

# 2D Ball2 radius =/= 1
b = Ball2([0., 0.], 2.)
# Test Dimension
@test dim(b) == 2
# Test Support Vector
d = [1., 0.]
@test σ(d, b) == [2., 0.]
d = [-1., 0.]
@test σ(d, b) == [-2., 0.]
d = [0., 1.]
@test σ(d, b) == [0., 2.]
d = [0., -1.]
@test σ(d, b) == [0., -2.]

# an_element function
b = Ball2([1., 2.], 2.)
@test an_element(b) ∈ b

# subset
b1 = Ball2([1., 2.], 2.)
b2 = Ball2([1., 2.], 0.)
b3 = Ball2([1.7, 2.7], 1.)
s = Singleton([1., 2.])
subset, point = ⊆(b1, s, true)
@test !⊆(b1, s) && !subset && point ∈ b1 && !(point ∈ s)
@test ⊆(b2, s) && ⊆(b2, s, true)[1]
@test !⊆(b1, BallInf([1., 2.], 1.))
@test ⊆(b2, BallInf([1., 2.], 2.))
subset, point = ⊆(b1, b3, true)
@test !⊆(b1, b3) && !subset && point ∈ b1 && !(point ∈ b3)
@test ⊆(b3, b1) && ⊆(b3, b1, true)[1]

# intersection
b1 = Ball2([0., 0.], 2.)
b2 = Ball2([2., 2.], 2.)
b3 = Ball2([4., 4.], 2.)
intersecting, point = is_intersection_empty(b1, b2, true)
@test is_intersection_empty(b1, b2) && intersecting && point ∈ b1 && point ∈ b2
@test !is_intersection_empty(b1, b3) && !is_intersection_empty(b1, b3, true)[1]
