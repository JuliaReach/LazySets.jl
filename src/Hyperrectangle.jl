import Base.LinAlg:norm

export Hyperrectangle, vertices_list, norm, radius, diameter, low, high

"""
    Hyperrectangle <: LazySet

Type that represents a hyperrectangle.

A [hyperrectangle](https://en.wikipedia.org/wiki/Hyperrectangle) is the Cartesian
product of one-dimensional intervals.

### Fields

- `center` -- center of the hyperrectangle as a real vector
- `radius` -- radius of the ball as a real vector, i.e., half of its width along
              each coordinate direction
"""
struct Hyperrectangle{N<:Real} <: LazySet
    center::Vector{N}
    radius::Vector{N}

    # default constructor
    Hyperrectangle{N}(center::Vector{N}, radius::Vector{N}) where {N<:Real} =
        (length(center) != length(radius)
            ? throw(DimensionMismatch)
            : new(center, radius))
end
# type-less convenience constructor
Hyperrectangle(center::Vector{N}, radius::Vector{N}) where {N<:Real} =
    Hyperrectangle{N}(center, radius)

"""
    Hyperrectangle(kwargs...)

Constructs a Hyperrectangle from keyword arguments.

### Input

Two combinations are allowed:

1. `center`, `radius` -- both vectors
2. `high`, `low`      -- both vectors (if both `center` and `radius` are also
                            defined, those are chosen instead)

### Examples

The following three constructions are equivalent:

```julia
julia> using LazySets
julia> c = ones(2);
julia> r = [0.1, 0.2];
julia> l = [0.9, 0.8];
julia> h = [1.1, 1.2];
julia> H1 = Hyperrectangle(c, r)
LazySets.Hyperrectangle([1.0, 1.0], [0.1, 0.2])
julia> H2 = Hyperrectangle(center=c, radius=r)
LazySets.Hyperrectangle([1.0, 1.0], [0.1, 0.2])
julia> H3 = Hyperrectangle(low=l, high=h)
LazySets.Hyperrectangle([1.0, 1.0], [0.1, 0.2])
```
"""
function Hyperrectangle(;kwargs...)
    dict = Dict{Symbol, Any}(kwargs)
    if length(dict) != 2
        # error below
    elseif haskey(dict, :center) && haskey(dict, :radius)
        return Hyperrectangle{eltype(dict[:center])}(dict[:center], dict[:radius])
    elseif haskey(dict, :high) && haskey(dict, :low)
        # compute center and radius from high and low vectors
        center = (dict[:high] .+ dict[:low]) ./ 2
        radius = abs.(dict[:high] .- center)
        return Hyperrectangle{eltype(center)}(center, radius)
    end
    throw(ArgumentError("Invalid arguments for Hyperrectangle: Use either " *
        "'center' and 'radius' or 'high' and 'low'."))
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
function σ(d::AbstractVector{<:Real}, H::Hyperrectangle)::AbstractVector{<:Real}
    return @. H.center + unit_step(d) * H.radius
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
function vertices_list(H::Hyperrectangle{N})::Vector{Vector{N}} where {N<:Real}
    return [H.center .+ si .* H.radius for si in IterTools.product([[1, -1] for i = 1:dim(H)]...)]
end

"""
    norm(H::Hyperrectangle, [p])

Return the norm of a Hyperrectangle. It is the norm of the enclosing ball (of
the given norm) of minimal volume.

### Input

- `H` -- hyperrectangle
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the norm.
"""
function norm(H::Hyperrectangle, p::Real=Inf)
    return maximum(map(x -> norm(x, p), vertices_list(H)))
end

"""
    radius(H::Hyperrectangle, [p])

Return the radius of a hyperrectangle. It is the radius of the enclosing ball
(of the given norm) of minimal volume with the same center.

### Input

- `H` -- hyperrectangle
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the radius.
"""
function radius(H::Hyperrectangle, p::Real=Inf)
    # the radius is the same for all corners of the hyperrectangle
    return norm(H.radius, p)
end

"""
    diameter(H::Hyperrectangle, [p])

Return the diameter of a hyperrectangle. It is the maximum distance between any
two elements of the set, or, equivalently, the diameter of the enclosing ball
(of the given norm) of minimal volume with the same center.

### Input

- `H` -- a hyperrectangle
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the diameter.
"""
function diameter(H::Hyperrectangle, p::Real=Inf)
    return radius(H, p) * 2
end

"""
    high(H::Hyperrectangle)

Return the higher coordinates of a hyperrectangle.

### Input

- `H` -- a hyperrectangle

### Output

A vector with the higher coordinates of the hyperrectangle, one entry per dimension.
"""
function high(H::Hyperrectangle)
    return H.center .+ H.radius
end

"""
    low(H::Hyperrectangle)

Return the lower coordinates of a hyperrectangle.

### Input

- `H` -- a hyperrectangle

### Output

A vector with the lower coordinates of the hyperrectangle, one entry per dimension.
"""
function low(H::Hyperrectangle)
    return H.center .- H.radius
end
