export Ball2

"""
    Ball2 <: LazySet

Type that represents a ball in the 2-norm.

### Fields

- `center` -- center of the ball as a real vector
- `radius` -- radius of the ball as a scalar (``≧ 0``)

### Examples

A five-dimensional ball in the 2-norm centered at the origin of radius 0.5:

```julia
julia> using LazySets
julia> B = Ball2(zeros(5), 0.5)
LazySets.Ball2([0.0, 0.0, 0.0, 0.0, 0.0], 0.5)
julia> dim(B)
5
```

We evaluate the support vector in a given direction:

```julia
julia> σ(ones(5), B)
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
        (radius < zero(N)
            ? throw(DomainError())
            : new(center, radius))
end
# type-less convenience constructor
Ball2(center::Vector{N}, radius::N) where {N<:Real} = Ball2{N}(center, radius)

"""
    dim(B)

Return the dimension of a Ball2.

### Input

- `B` -- a ball in the 2-norm

### Output

The ambient dimension of the ball.
"""
function dim(B::Ball2)::Int64
    length(B.center)
end

"""
    σ(d, B)

Return the support vector of a Ball2 in a given direction.

### Input

- `d` -- a direction
- `B` -- a ball in the 2-norm

### Output

The support vector in the given direction.

### Notes

If the given direction has norm zero, the origin is returned.
"""
function σ(d::AbstractVector{<:Real}, B::Ball2)::Vector{<:Real}
    dnorm = norm(d)
    if dnorm > 0
        return B.center .+ d .* (B.radius / dnorm)
    else
        return zeros(eltype(d), length(d))
    end
end

