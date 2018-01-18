# constructor from empty set
@test SymmetricIntervalHull(∅) == ∅

# normal constructor
h = SymmetricIntervalHull(Ball2([2., 3.], 4.))

# dimension
@test dim(h) == 2

# support vector
d = [1., 0.]
@test σ(d, h)[1] == 6. && σ(-d, h)[1] == -6.
d = [0., 1.]
@test σ(d, h)[2] == 7. && σ(-d, h)[2] == -7.

# an_element function
@test an_element(h) == [0., 0.]
