# normal constructor
normal = ones(3)
h = HalfSpace(normal, 5.)

# dimension
@test dim(h) == 3

# support vector and membership function
function test_svec(h, d)
    @test σ(d, h) ∈ h
    @test σ(2. * d, h) ∈ h
    d = [1., 0., 0.]
    @test_throws ErrorException σ(d, h)
    d = zeros(3)
    @test σ(d, h) ∈ h
end
# tests 1
normal = ones(3)
d = ones(3)
test_svec(HalfSpace(normal, 5.), d)
# tests 2
normal = zeros(3); normal[3] = 1.
d = zeros(3); d[3] = 1.
test_svec(HalfSpace(normal, 5.), d)

# an_element function and membership function
@test an_element(h) ∈ h
