# normal constructor
normal = ones(3)
hs = HalfSpace(normal, 5.)

# dimension
@test dim(hs) == 3

# support vector and membership function
function test_svec(hs, d)
    @test σ(d, hs) ∈ hs
    @test σ(2. * d, hs) ∈ hs
    d2 = [1., 0., 0.]
    @test_throws ErrorException σ(d2, hs)
    d2 = zeros(3)
    @test σ(d2, hs) ∈ hs
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
@test an_element(hs) ∈ hs
