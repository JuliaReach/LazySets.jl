# 1D ellipsoid
E = Ellipsoid([0.], diagm(1.))
# Test Dimension
@test dim(E) == 1
# Test Support Vector
d = [1.]
@test σ(d, E) == [1.]
d = [-1.]
@test σ(d, E) == [-1.]
# test constructor
E = Ellipsoid(diagm(1.))
@test E.center == [0.]

# 2D Ellipsoid
E = Ellipsoid([0., 0.], diagm([1., 1.]))
# Test Dimension
@test dim(E) == 2
# Test Support Vector
d = [1., 0.]
@test σ(d, E) == [1., 0.]
d = [-1., 0.]
@test σ(d, E) == [-1., 0.]
d = [0., 1.]
@test σ(d, E) == [0., 1.]
d = [0., -1.]
@test σ(d, E) == [0., -1.]

# 2D Ellipsoid not 0-centered
E = Ellipsoid([1., 2.], diagm([1., 1.]))
# Test Dimension
@test dim(E) == 2
# Test Support Vector
d = [1., 0.]
@test σ(d, E) == [2., 2.]
d = [-1., 0.]
@test σ(d, E) == [0., 2.]
d = [0., 1.]
@test σ(d, E) == [1., 3.]
d = [0., -1.]
@test σ(d, E) == [1., 1.]

# another shape matrix
E = Ellipsoid([1., 2.], diagm([.5, 2.]))
# Test Support Vector
d = [1., 0.]
@test σ(d, E) ≈ [1+sqrt(.5), 2.]
d = [-1., 0.]
@test σ(d, E) ≈ [1-sqrt(.5), 2.]
d = [0., 1.]
@test σ(d, E) ≈ [1., 2. + sqrt(2)]
d = [0., -1.]
@test σ(d, E) ≈ [1., 2. - sqrt(2)]

# an_element and set membership functions
E = Ellipsoid([1., 2.], 2*eye(2))
@test an_element(E) ∈ E
