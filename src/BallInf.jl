import Base.LinAlg:norm

export BallInf, vertices_list, norm, radius, diameter

"""
    BallInf <: LazySet

Type that represents a ball in the infinity norm.

### Fields

- `center` -- center of the ball as a real vector
- `radius` -- radius of the ball as a scalar (``≧ 0``)

### Examples

We create the two-dimensional unit ball, and compute its support function
along the direction ``(1, 1)``:

```julia
julia> B = BallInf(zeros(2), 0.1)
LazySets.BallInf([0.0, 0.0], 0.1)

julia> dim(B)
2

julia> ρ([1., 1.], B)
0.2
```
"""
struct BallInf{N<:Real} <: LazySet
    center::Vector{N}
    radius::N

    # default constructor with domain constraint for radius
    BallInf{N}(center, radius) where N =
        (radius < zero(N)
            ? throw(DomainError())
            : new(center, radius))
end
# type-less convenience constructor
BallInf(center::Vector{N}, radius::N) where {N<:Real} = BallInf{N}(center, radius)

"""
    dim(B)

Return the dimension of a BallInf.

### Input

- `B` -- a ball in the infinity norm

### Output

The ambient dimension of the ball.
"""
function dim(B::BallInf)::Int64
    return length(B.center)
end

"""
    σ(d, B)

Return the support vector of an infinity-norm ball in a given direction.

### Input

- `d` -- direction
- `B` -- unit ball in the infinity norm

### Algorithm

This code is a vectorized version of

```julia
[(d[i] >= 0) ? B.center[i] + B.radius : B.center[i] - B.radius for i in 1:length(d)]
```

Notice that we cannot use `B.center + sign.(d) * B.radius`, since the built-in `sign`
function is such that `sign(0) = 0`, instead of 1. For this reason, we use the
custom `unit_step` function, that allows to do: `B.center + unit_step.(d) * B.radius`
(the dot operator performs broadcasting, to accept vector-valued entries).
"""
function σ(d::AbstractVector{<:Real}, B::BallInf)::Vector{<:Real}
    return B.center .+ unit_step.(d) .* B.radius
end

"""
    vertices_list(B::BallInf)

Return the list of vertices of a ball in the infinity norm.

### Input

- `B` -- a ball in the infinity norm

### Output

The list of vertices as an array of floating-point vectors.

### Notes

For high-dimensions, it is preferable to develop a `vertex_iterator` approach.
"""
function vertices_list(B::BallInf{N})::Vector{Vector{N}} where {N<:Real}
    return [B.center .+ si .* B.radius for si in IterTools.product([[1, -1] for i = 1:dim(B)]...)]
end

"""
    norm(B::BallInf, [p])

Return the norm of a `BallInf`. It is the norm of the enclosing ball (of
the given norm) of minimal volume.

### Input

- `B` -- ball in the infinity norm
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the norm.
"""
function norm(B::BallInf, p::Real=Inf)::Real
    return maximum(map(x -> norm(x, p), vertices_list(B)))
end

"""
    radius(B::BallInf, [p])

Return the radius of a ball in the infinity norm. It is the radius of the
enclosing ball (of the given norm) of minimal volume with the same center.

### Input

- `B` -- a ball in the infinity norm
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the radius.
"""
function radius(B::BallInf, p::Real=Inf)::Real
    if p == Inf
        return B.radius
    else
        return norm(fill(B.radius, dim(B)), p)
    end
end

"""
    diameter(B::BallInf, [p])

Return the diameter of a ball in the infinity norm. It is the maximum distance
between any two elements of the set, or, equivalently, the diameter of the
enclosing ball (of the given norm) of minimal volume with the same center.

### Input

- `B` -- a ball in the infinity norm
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the diameter.
"""
function diameter(B::BallInf, p::Real=Inf)::Real
    return radius(B, p) * 2
end
