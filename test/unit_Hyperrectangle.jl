# 1D Hyperrectangle
h = Hyperrectangle([0.], [1.])
# Test Dimension
@test dim(h) == 1
# Test Support Vector
d = [1.]
@test σ(d, h) == [1.]
d = [-1.]
@test σ(d, h) == [-1.]

# 2D Hyperrectangle
h = Hyperrectangle([0., 0.], [1., 1.])
# Test Dimension
@test dim(h) == 2
# Test Support Vector
d = [1., 1.]
@test σ(d, h) == [1., 1.]
d = [-1., 1.]
@test σ(d, h) == [-1., 1.]
d = [-1., -1.]
@test σ(d, h) == [-1., -1.]
d = [1., -1.]
@test σ(d, h) == [1., -1.]

# 2D Hyperrectangle not 0-centered
h = Hyperrectangle([1., 2.], [1., 1.])
# Test Dimension
@test dim(h) == 2
# Test Support Vector
d = [1., 1.]
@test σ(d, h) == [2., 3.]
d = [-1., 1.]
@test σ(d, h) == [0., 3.]
d = [-1., -1.]
@test σ(d, h) == [0., 1.]
d = [0., -1.]
@test σ(d, h) == [2., 1.]

# 2D Hyperrectangle not same radius in each direction
h = Hyperrectangle([0., 0.], [1., 2.])
# Test Dimension
@test dim(h) == 2
# Test Support Vector
d = [1., 1.]
@test σ(d, h) == [1., 2.]
d = [-1., 1.]
@test σ(d, h) == [-1., 2.]
d = [-1., -1.]
@test σ(d, h) == [-1., -2.]
d = [1., -1.]
@test σ(d, h) == [1., -2.]

function isin(e, list)
    for x in list
        if x == e
            return true
        end
    end
    return false
end

# 2D Hyperrectangle not centered, not same radius, for vertex representation,
# radius, and diameter
h = Hyperrectangle([3., 2.], [2., 1.])
vl = vertices_list(h)
# Test Vertices
@test length(vl) == 4
@test isin([1., 1.], vl) == true
@test isin([1., 3.], vl) == true
@test isin([5., 1.], vl) == true
@test isin([5., 3.], vl) == true
# norm
@test norm(h) == norm([5., 3.], Inf)
# radius
@test radius(h) == norm([2., 1.], Inf)
# diameter
@test diameter(h) == norm([5., 3.] - [1., 1.], Inf)

# alternative constructors
c = ones(2)
r = [0.1, 0.2]
l = [0.9, 0.8]
h = [1.1, 1.2]
H1 = Hyperrectangle(c, r)
H2 = Hyperrectangle(center=c, radius=r)
H3 = Hyperrectangle(low=l, high=h)
@test H1.center == H2.center
@test H2.center ≈ H3.center
@test H1.radius == H2.radius
@test H2.radius ≈ H3.radius

# Test low and high methods for a hyperrectangle
H = Hyperrectangle(center=[-2.1, 5.6, 0.9], radius=fill(0.5, 3))
@test low(H) == [-2.6, 5.1, 0.4]
@test high(H) == [-1.6, 6.1, 1.4]

# membership
H = Hyperrectangle([1.0, 1.0], [2.0, 3.0])
@test !∈([-1.1, 4.1], H)
@test ∈([-1.0, 4.0], H)

# an_element function
H = Hyperrectangle([1.0, 2.0], [3.0, 4.0])
@test an_element(H) ∈ H

# subset
H1 = Hyperrectangle([1.5, 1.5], [0.5, 0.5])
H2 = Hyperrectangle([2.0, 2.5], [0.5, 0.5])
H3 = Hyperrectangle([2.0, 2.0], [2.0, 3.0])
B1 = BallInf([2.0, 2.5], 0.5)
B2 = BallInf([2.0, 2.0], 1.0)
@test !⊆(H1, H2) && ⊆(H1, H3) && ⊆(H2, H3)
@test ⊆(H2, B1) && ⊆(B1, H2)
@test ⊆(B1, B2) && !⊆(B2, B1)
