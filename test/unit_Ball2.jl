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
@test !⊆(Ball2([1., 2.], 2.), Singleton([1., 2.])) && ⊆(Ball2([1., 2.], 0.), Singleton([1., 2.]))
