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
