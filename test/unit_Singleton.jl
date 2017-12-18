# 1D singleton
s = Singleton([1.])
d = [1.]
@test σ(d, s) == [1.]
d = [-1.]
@test σ(d, s) == [1.]
d = [0.]
@test σ(d, s) == [1.]

# 2D singleton
s = Singleton([1., 2.])
d = [1., 0.]
@test σ(d, s) == [1., 2.]
d = [-1., .5]
@test σ(d, s) == [1., 2.]
d = [0., 0.]
@test σ(d, s) == [1., 2.]

# membership
S = Singleton([1., 1.])
!∈([0.9, 1.1], S)
∈([1.0, 1.0], S)

# an_element function
@test an_element(S) ∈ S

# subset
@test ⊆(Singleton([0., 1.]), VPolygon([[0.,0.],[0., 2.]]))
@test !⊆(Singleton([0., 3.]), VPolygon([[0.,0.],[0., 2.], [2., 0.]]))
