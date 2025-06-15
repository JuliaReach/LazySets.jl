"""
    Ball1{N, VN<:AbstractVector{N}} <: AbstractCentrallySymmetricPolytope{N}

Type that represents a ball in the 1-norm (also known as the Manhattan norm).
The ball is also known as a
[cross-polytope](https://en.wikipedia.org/wiki/Cross-polytope).

It is defined as the set

```math
\\mathcal{B}_1^n(c, r) = \\{ x ∈ ℝ^n : ∑_{i=1}^n |c_i - x_i| ≤ r \\},
```
where ``c ∈ ℝ^n`` is its center and ``r ∈ ℝ_+`` its radius.

### Fields

- `center` -- center of the ball as a real vector
- `radius` -- radius of the ball as a scalar (``≥ 0``)

### Examples

The unit ball in the 1-norm in the plane:

```jldoctest ball1_constructor
julia> B = Ball1(zeros(2), 1.0)
Ball1{Float64, Vector{Float64}}([0.0, 0.0], 1.0)
julia> dim(B)
2
```

We evaluate the support vector in the North direction:

```jldoctest ball1_constructor
julia> σ([0.0, 1.0], B)
2-element Vector{Float64}:
 0.0
 1.0
```
"""
struct Ball1{N,VN<:AbstractVector{N}} <: AbstractCentrallySymmetricPolytope{N}
    center::VN
    radius::N

    # default constructor with domain constraint for radius
    function Ball1(center::VN, radius::N) where {N,VN<:AbstractVector{N}}
        @assert radius >= zero(N) "the radius must be nonnegative but is $radius"
        return new{N,VN}(center, radius)
    end
end
