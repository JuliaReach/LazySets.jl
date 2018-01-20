# normal constructor
normal = ones(3)
hp = Hyperplane(normal, 5.)

# dimension
@test dim(hp) == 3

# support vector and membership function
function test_svec(hp, d)
    @test σ(d, hp) ∈ hp
    @test σ(2. * d, hp) ∈ hp
    d2 = [1., 0., 0.]
    @test_throws ErrorException σ(d2, hp)
    d2 = zeros(3)
    @test σ(d2, hp) ∈ hp
end
# tests 1
normal = ones(3)
d = ones(3)
test_svec(Hyperplane(normal, 5.), d)
# tests 2
normal = zeros(3); normal[3] = 1.
d = zeros(3); d[3] = 1.
test_svec(Hyperplane(normal, 5.), d)

# an_element function and membership function
@test an_element(hp) ∈ hp
