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

# 2D Hyperrectangle not centered, not same radius, for vertex representation, radius,
# and diameter
h = Hyperrectangle([3., 2.], [2., 1.])
vl = vertices_list(h)
# Test Vertices
@test length(vl) == 4
@test isin([1., 1.], vl) == true
@test isin([1., 3.], vl) == true
@test isin([5., 1.], vl) == true
@test isin([5., 3.], vl) == true
# Test Radius
@test radius(h) == norm([5., 3.])
# Test Diameter
@test diameter(h) == norm([5., 3.] - [1., 1.])
