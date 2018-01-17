B = BallInf(ones(2), 3.)
H = Hyperrectangle(ones(2), ones(2))
I = Intersection(B, H)

# dim
@test dim(I) == 2

# support vector (error)
@test_throws ErrorException σ(ones(2), I)

# emptiness of intersection
@test !isempty(I)
