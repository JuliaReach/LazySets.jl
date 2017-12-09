export Ball1

"""
    Ball1 <: LazySet

Type that represents a ball in the 1-norm, also known as Manhattan or Taxicab
norm.

It is defined as the set

```math
\\mathcal{B}_1^n(c, r) = \\{ x ∈ \\mathbb{R}^n : ∑_{i=1}^n |x_i| ≤ r \\},
```
where ``c ∈ \\mathbb{R}^n`` is its center and ``r ∈ \\mathbb{R}_+`` its radius.

### Fields

- `center` -- center of the ball as a real vector
- `radius` -- radius of the ball as a scalar (``≥ 0``)

### Examples

Unit ball in the 1-norm in the plane:

```jldoctest ball1_constructor
julia> B = Ball1(zeros(2), 1.)
LazySets.Ball1{Float64}([0.0, 0.0], 1.0)
julia> dim(B)
2
```

We evaluate the support vector in the East direction:

```jldoctest ball1_constructor
julia> σ([0.,1], B)
1-element Array{Float64,1}:
 1.0
 0.0
```
"""
struct Ball1{N<:Real} <: LazySet
    center::Vector{N}
    radius::N

    # default constructor with domain constraint for radius
    Ball1{N}(center, radius) where N =
        radius < zero(N) ? throw(DomainError()) : new(center, radius)
end
# type-less convenience constructor
Ball1(center::Vector{N}, radius::N) where {N<:Real} = Ball1{N}(center, radius)

"""
    dim(B::Ball1)::Int

Return the dimension of a `Ball1`.

### Input

- `B` -- a ball in the 1-norm

### Output

The ambient dimension of the ball.
"""
function dim(B::Ball1)::Int
    return length(B.center)
end

"""
    σ(d::AbstractVector{N}, B::Ball1)::AbstractVector{N} where {N<:AbstractFloat}

Return the support vector of a `Ball1` in a given direction.

### Input

- `d` -- a direction
- `B` -- a ball in the p-norm

### Output

Support vector in the given direction.
"""
function σ(d::AbstractVector{N},
           B::Ball1)::AbstractVector{N} where {N<:AbstractFloat}
    res = copy(B.center)
    imax = indmax(abs.(d)) 
    res[imax] = sign(d[imax]) * B.radius
    return res
end
