export Ball2

"""
    Ball2 <: LazySet

Type that represents a ball in the 2-norm.

It is defined as the set

```math
\\mathcal{B}_2^n(c, r) = \\{ x ∈ \\mathbb{R}^n : ‖ x - c ‖_2 ≤ r \\},
```
where ``c ∈ \\mathbb{R}^n`` is its center and ``r ∈ \\mathbb{R}_+`` its radius.
Here ``‖ ⋅ ‖_2`` denotes the Euclidean norm (also known as 2-norm), defined as
``‖ x ‖_2 = \\left( \\sum\\limits_{i=1}^n |x_i|^2 \\right)^{1/2}`` for any
``x ∈ \\mathbb{R}^n``.

### Fields

- `center` -- center of the ball as a real vector
- `radius` -- radius of the ball as a real scalar (``≥ 0``)

### Examples

A five-dimensional ball in the 2-norm centered at the origin of radius 0.5:

```julia
julia> using LazySets
julia> B = Ball2(zeros(5), 0.5)
LazySets.Ball2([0.0, 0.0, 0.0, 0.0, 0.0], 0.5)
julia> dim(B)
5
```

Evaluate its support vector in the direction ``[1,2,3,4,5]``:

```julia
julia> σ([1.,2.,3.,4.,5.], B)
5-element Array{Float64,1}:
 0.06742
 0.13484
 0.20226
 0.26968
 0.3371
```
"""
struct Ball2{N<:Real} <: LazySet
    center::Vector{N}
    radius::N

    # default constructor with domain constraint for radius
    Ball2{N}(center, radius) where N =
        (radius < zero(N) ? throw(DomainError()) : new(center, radius))
end
# type-less convenience constructor
Ball2(center::Vector{N}, radius::N) where {N<:Real} = Ball2{N}(center, radius)

"""
    dim(B)

Return the dimension of a ball in the 2-norm.

### Input

- `B` -- a ball in the 2-norm

### Output

The ambient dimension of the ball.
"""
dim(B::Ball2) = length(B.center)

"""
    σ(d, B)

Return the support vector of a ball in the 2-norm in a given direction.

### Input

- `d` -- a direction
- `B` -- a ball in the 2-norm

### Output

The support vector in the given direction. If the given direction has norm zero,
the origin is returned.

### Notes

This function requires computing the 2-norm of the input direction, which is
performed in the given precision of the numeric datatype of both the direction
and the set.
Exact inputs are not supported.
"""
function σ(d::AbstractVector{N},
           B::Ball2)::AbstractVector{<:AbstractFloat} where {N<:AbstractFloat}
    dnorm = norm(d, 2)
    if dnorm <= zero(N)
        return zeros(eltype(d), length(d))
    else
        return @. B.center + d * (B.radius / dnorm)
    end
end
