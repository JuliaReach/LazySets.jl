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
s1 = Singleton([0., 1.])
s2 = Singleton([0., 3.])
p1 = VPolygon([[0.,0.],[0., 2.]])
p2 = VPolygon([[0.,0.],[0., 2.], [2., 0.]])
@test ⊆(s1, p1) && ⊆(s1, p1, true)[1]
subset, point = ⊆(s2, p2, true)
@test !⊆(s2, p2) && !subset && point ∈ s2 && !(point ∈ p2)
@test ⊆(s1, s1) && ⊆(s1, s1, true)[1]
subset, point = ⊆(s1, s2, true)
@test !⊆(s1, s2) && !subset && point ∈ s1 && !(point ∈ s2)

# intersection
S1 = Singleton([1.0, 1.0])
S2 = Singleton([0.0, 0.0])
S3 = ZeroSet(2)
@test !∩(S1, S2) && !∩(S1, S2, true)[1]
intersection, point = ∩(S2, S3, true)
@test ∩(S2, S3) && intersection && point ∈ S2 && point ∈ S3
