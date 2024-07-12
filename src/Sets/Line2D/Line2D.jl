"""
    Line2D{N, VN<:AbstractVector{N}} <: AbstractPolyhedron{N}

Type that represents a line in 2D of the form ``a⋅x = b`` (i.e., a special case
of a `Hyperplane`).

### Fields

- `a` -- normal direction (non-zero)
- `b` -- constraint

### Examples

The line ``y = -x + 1``:

```jldoctest
julia> Line2D([1., 1.], 1.)
Line2D{Float64, Vector{Float64}}([1.0, 1.0], 1.0)
```

The alternative constructor takes two 2D points (`AbstractVector`s) `p` and `q`
and creates a canonical line from `p` to `q`. See the algorithm section below
for details.

```jldoctest
julia> Line2D([1., 1.], [2., 2])
Line2D{Float64, Vector{Float64}}([-1.0, 1.0], 0.0)
```

### Algorithm

Given two points ``p = (x₁, y₁)`` and ``q = (x₂, y₂)``, the line that passes
through these points is

```math
ℓ:~~y - y₁ = \\dfrac{(y₂ - y₁)}{(x₂ - x₁)} ⋅ (x-x₁).
```
The particular case ``x₂ = x₁`` defines a line parallel to the ``y``-axis
(vertical line).
"""
struct Line2D{N,VN<:AbstractVector{N}} <: AbstractPolyhedron{N}
    a::VN
    b::N

    # default constructor with length constraint
    function Line2D(a::VN, b::N) where {N,VN<:AbstractVector{N}}
        @assert length(a) == 2 "a Line2D must be two-dimensional"
        @assert !iszero(a) "a line needs a non-zero normal vector"
        return new{N,VN}(a, b)
    end
end

function Line2D(p::AbstractVector, q::AbstractVector)
    @assert length(p) == length(q) == 2 "a Line2D must be two-dimensional"

    N = promote_type(eltype(p), eltype(q))
    x₁, y₁ = @inbounds p[1], p[2]
    x₂, y₂ = @inbounds q[1], q[2]

    if x₁ == x₂  # line is vertical
        @assert y₁ != y₂ "a line needs two distinct points"
        a = [one(N), zero(N)]
        b = x₁
        return Line2D(a, b)
    end

    k = (y₁ - y₂) / (x₂ - x₁)
    a = [k, one(N)]
    b = y₁ + k * x₁
    return Line2D(a, b)
end
