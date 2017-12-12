import Base.LinAlg:norm,
       Base.Iterators.repeated

export BallInf,
       vertices_list,
       norm,
       radius,
       diameter

"""
    BallInf{N<:Real} <: LazySet

Type that represents a ball in the infinity norm.

### Fields

- `center` -- center of the ball as a real vector
- `radius` -- radius of the ball as a real scalar (``≥ 0``)

### Notes

Mathematically, a ball in the infinity norm is defined as the set

```math
\\mathcal{B}_∞^n(c, r) = \\{ x ∈ \\mathbb{R}^n : ‖ x - c ‖_∞ ≤ r \\},
```
where ``c ∈ \\mathbb{R}^n`` is its center and ``r ∈ \\mathbb{R}_+`` its radius.
Here ``‖ ⋅ ‖_∞`` denotes the infinity norm, defined as
``‖ x ‖_∞ = \\max\\limits_{i=1,…,n} \\vert x_i \\vert`` for any
``x ∈ \\mathbb{R}^n``.

### Examples

Create the two-dimensional unit ball and compute its support function along the
positive ``x=y`` direction:

```jldoctest
julia> B = BallInf(zeros(2), 1.0)
LazySets.BallInf{Float64}([0.0, 0.0], 1.0)
julia> dim(B)
2
julia> ρ([1., 1.], B)
2.0
```
"""
struct BallInf{N<:Real} <: LazySet
    center::Vector{N}
    radius::N

    # default constructor with domain constraint for radius
    BallInf{N}(center, radius) where N =
        radius < zero(N) ? throw(DomainError()) : new(center, radius)
end
# type-less convenience constructor
BallInf(center::Vector{N}, radius::N) where {N<:Real} =
    BallInf{N}(center, radius)

"""
    dim(B::BallInf)::Int

Return the dimension of a ball in the infinity norm.

### Input

- `B` -- ball in the infinity norm

### Output

The ambient dimension of the ball.
"""
function dim(B::BallInf)::Int
    return length(B.center)
end

"""
    σ(d::AbstractVector{<:Real}, B::BallInf)::AbstractVector{<:Real}

Return the support vector of an infinity norm ball in a given direction.

### Input

- `d` -- direction
- `B` -- ball in the infinity norm

### Output

The support vector in the given direction.
If the direction has norm zero, the vertex with biggest values is returned.
"""
function σ(d::AbstractVector{<:Real}, B::BallInf)::AbstractVector{<:Real}
    return @. B.center + sign_cadlag(d) * B.radius
end

"""
    vertices_list(B::BallInf{N})::Vector{Vector{N}} where {N<:Real}

Return the list of vertices of a ball in the infinity norm.

### Input

- `B` -- ball in the infinity norm

### Output

A list of vertices.

### Notes

For high dimensions, it is preferable to develop a `vertex_iterator` approach.

### Examples

```jldoctest
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
function vertices_list(B::BallInf{N})::Vector{Vector{N}} where {N<:Real}
    return [B.center .+ si .* B.radius
        for si in IterTools.product([[1, -1] for i = 1:dim(B)]...)]
end

"""
    norm(B::BallInf, [p]::Real=Inf)::Real

Return the norm of a ball in the infinity norm.

### Input

- `B` -- ball in the infinity norm
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the norm.

### Notes

The norm of an infinity ball is defined as the norm of the enclosing ball, of
the given ``p``-norm, of minimal volume.
"""
function norm(B::BallInf, p::Real=Inf)::Real
    return maximum(map(x -> norm(x, p), vertices_list(B)))
end

"""
    radius(B::BallInf, [p]::Real=Inf)::Real

Return the radius of a ball in the infinity norm.

### Input

- `B` -- ball in the infinity norm
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the radius.

### Notes

The radius is defined as the radius of the enclosing ball of the given
``p``-norm of minimal volume with the same center.
"""
function radius(B::BallInf, p::Real=Inf)::Real
    return (p == Inf) ? B.radius : norm(fill(B.radius, dim(B)), p)
end

"""
    diameter(B::BallInf, [p]::Real=Inf)::Real

Return the diameter of a ball in the infinity norm.

### Input

- `B` -- ball in the infinity norm
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the diameter.

### Notes

The diameter is defined as the maximum distance in the given ``p``-norm between
any two elements of the set.
Equivalently, it is the diameter of the enclosing ball of the given ``p``-norm
of minimal volume with the same center.
"""
function diameter(B::BallInf, p::Real=Inf)::Real
    return radius(B, p) * 2
end
