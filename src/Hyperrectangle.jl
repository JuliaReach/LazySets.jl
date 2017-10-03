"""
    Hyperrectangle <: LazySet

Type that represents a Hyperrectangle.

FIELDS:

- ``c`` -- the center
- ``r`` -- the radius, i.e. the width of the rectangle in each direction
"""
struct Hyperrectangle <: LazySet
    center::Vector{Float64}
    radius::Vector{Float64}
    Hyperrectangle(center, radius) = length(center) != length(radius) ? throw(DimensionMismatch) : new(center, radius)
end

"""
    dim(H)

Return the dimension of a Hyperrectangle.

INPUT:

- ``H`` -- a hyperrectangle

OUTPUT:

The ambient dimension of the hyperrectangle.
"""
function dim(H::Hyperrectangle)::Int64
    length(H.center)
end

"""
    σ(d, H)

Return the support vector of a Hyperrectangle in a given direction.

See also: ``BallInf``.
"""
function σ(d::Union{Vector{Float64}, SparseVector{Float64,Int64}}, H::Hyperrectangle)::Vector{Float64}
    return H.center .+ unit_step.(d) .* H.radius
end

using IterTools

"""
    vertices_list(H)

Return the vertices of a hyperrectangle.

INPUT:

- ``H`` -- a hyperrectangle

OUTPUT:

The list of vertices as an array of floating-point vectors.

NOTES:

For high-dimensions, it is preferable to develop a ``vertex_iterator`` approach.
"""
function vertices_list(H::Hyperrectangle)::Array{Vector{Float64}, 1}
    return [H.center .+ si .* H.radius for si in IterTools.product([[1, -1] for i = 1:dim(H)]...)]
end

"""
    radius(H)

Return the radius of a Hyperrectangle. It is the radius of the enclosing
hypercube of minimal volume.

INPUT:

- ``H`` -- a hyperrectangle

OUTPUT:

A real number representing its radius.
"""
function radius(H::Hyperrectangle)::Float64
    return maximum(map(norm, vertices_list(H)))
end

"""
    diameter(H)

Return the diameter of a hyperrectangle. It the maximum norm (measured the
infinity norm) of any element of the set.

INPUT:

- ``H`` -- a hyperrectangle

OUTPUT:

The diameter of the hyperrectangle.
"""
function diameter(H::Hyperrectangle)::Float64
    return 2. * radius(Hyperrectangle(zeros(dim(H)), H.radius))
end

export Hyperrectangle, vertices_list, radius, diameter
