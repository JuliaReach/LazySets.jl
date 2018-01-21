B = BallInf(ones(2), 3.)
H = Hyperrectangle(ones(2), ones(2))
I = Intersection(B, H)

# dim
@test dim(I) == 2

# support vector (error)
@test_throws ErrorException σ(ones(2), I)

# membership
@test ∈(ones(2), I) && !∈([5., 5.], I)

# emptiness of intersection
@test !isempty(I)
