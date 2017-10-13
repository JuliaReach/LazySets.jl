using IterTools

export Hyperrectangle, vertices_list, radius, diameter

"""
    Hyperrectangle <: LazySet

Type that represents a Hyperrectangle.

A [hyperrectangle](https://en.wikipedia.org/wiki/Hyperrectangle) is the Cartesian
product of one-dimensional intervals.

### Fields

- `center` -- center of the hyperrectangle as a real vector
- `radius` -- radius of the ball as a real vector, i.e. its width along
              each coordinate direction
"""
struct Hyperrectangle <: LazySet
    center::Vector{Float64}
    radius::Vector{Float64}
    Hyperrectangle(center, radius) = length(center) != length(radius) ?
                                throw(DimensionMismatch) : new(center, radius)
end

"""
    dim(H)

Return the dimension of a Hyperrectangle.

### Input

- `H` -- a hyperrectangle

### Output

The ambient dimension of the hyperrectangle as an integer.
"""
function dim(H::Hyperrectangle)::Int64
    length(H.center)
end

"""
    σ(d, H)

Return the support vector of a Hyperrectangle in a given direction.
"""
function σ(d::Union{Vector{Float64}, SparseVector{Float64,Int64}},
           H::Hyperrectangle)::Vector{Float64}
    return H.center .+ unit_step.(d) .* H.radius
end

"""
    vertices_list(H::Hyperrectangle)

Return the vertices of a hyperrectangle.

### Input

- `H` -- a hyperrectangle

### Output

The list of vertices as an array of floating-point vectors.

### Notes

For high-dimensions, it is preferable to develop a `vertex_iterator` approach.
"""
function vertices_list(H::Hyperrectangle)::Array{Vector{Float64}, 1}
    return [H.center .+ si .* H.radius for si in IterTools.product([[1, -1] for i = 1:dim(H)]...)]
end

"""
    radius(H::Hyperrectangle)

Return the radius of a Hyperrectangle. It is the radius of the enclosing
hypercube of minimal volume.

### Input

- `H` -- a hyperrectangle

### Output

A real number representing its radius.
"""
function radius(H::Hyperrectangle)::Float64
    return maximum(map(norm, vertices_list(H)))
end

"""
    diameter(H::Hyperrectangle)

Return the diameter of a hyperrectangle. It the maximum norm (measured the
infinity norm) of any element of the set.

### Input

- `H` -- a hyperrectangle

### Output

The diameter of the hyperrectangle.
"""
function diameter(H::Hyperrectangle)::Float64
    return 2. * radius(Hyperrectangle(zeros(dim(H)), H.radius))
end
