import Base.LinAlg:norm

export BallInf, vertices_list, norm, radius, diameter

"""
    BallInf <: LazySet

Type that represents a ball in the infinity norm.

It is defined as the set

```math
\\mathcal{B}_∞^n(c, r) = \\{ x ∈ \\mathbb{R}^n : ‖ x - c ‖_∞ ≦ r \\},
```
where ``c ∈ \\mathbb{R}^n`` is its center and ``r ∈ \\mathbb{R}_+`` its radius.
Here ``‖ ⋅ ‖_∞`` denotes the infinity norm, defined as
``‖ x ‖_∞ = \\max\\limits_{i=1,…,n} \\vert x_i \\vert`` for any
``x ∈ \\mathbb{R}^n``.

### Fields

- `center` -- center of the ball as a real vector
- `radius` -- radius of the ball as a scalar (``≧ 0``)

### Examples

We create the two-dimensional unit ball, and compute its support function
along the positive ``x=y`` direction:

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
    dim(B::BallInf)

Return the dimension of a ball in the infinity norm.

### Input

- `B` -- a ball in the infinity norm

### Output

The ambient dimension of the ball.
"""
dim(B::BallInf) = length(B.center)

"""
    σ(d::AbstractVector{<:Real}, B::BallInf)

Return the support vector of an infinity norm ball in a given direction.

### Input

- `d` -- direction
- `B` -- unit ball in the infinity norm
"""
function σ(d::AbstractVector{<:Real}, B::BallInf)::AbstractVector{<:Real}
    #=
    We cannot use `B.center + sign.(d) * B.radius`, since the built-in `sign`
    function is such that `sign(0) = 0`, instead of 1 as needed.
    The custom `unit_step` function allows to do perform broadcasting through
    the dot operator, thus accepting vector-valued entries.
    =#
    return @. B.center + unit_step(d) * B.radius
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

### Examples

```julia
julia> B = BallInf(zeros(2), 0.1)
LazySets.BallInf{Float64}([0.0, 0.0], 0.1)
julia> vertices_list(B)
4-element Array{Array{Float64,1},1}:
 [0.1, 0.1]
 [-0.1, 0.1]
 [0.1, -0.1]
 [-0.1, -0.1]
```
"""
function vertices_list(B::BallInf)
    return [B.center .+ si .* B.radius for si in IterTools.product([[1, -1] for i = 1:dim(B)]...)]
end

"""
    norm(B::BallInf, [p])

Return the norm of a ball in the infinity norm.

It is defined as the norm of the enclosing ball, of the given
``p``-norm, of minimal volume.

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

Return the radius of a ball in the infinity norm.

It is defined as the radius of the enclosing ball of the given ``p``-norm of
minimal volume with the same center.

### Input

- `B` -- a ball in the infinity norm
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the radius.
"""
radius(B::BallInf, p::Real=Inf) = (p == Inf) ? B.radius : norm(fill(B.radius, dim(B)), p)

"""
    diameter(B::BallInf, [p])

Return the diameter of a ball in the infinity norm.

It corresponds to the maximum distance between any two elements of the set.
Equivalently, it is the diameter of the enclosing ball of the given ``p``-norm
of minimal volume with the same center.

### Input

- `B` -- a ball in the infinity norm
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the diameter.
"""
diameter(B::BallInf, p::Real=Inf) = 2 * radius(B, p)
