export Ball2

"""
    Ball2 <: LazySet

Type that represents a ball in the 2-norm.

This set is defined by a center and a radius,

    ``B_2(c, r) = \\{ x ∈ \\mathbb{R}^n : ‖x - c‖_2 ≦ r \\}.``

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
struct Ball2 <: LazySet
    center::Vector{Float64}
    radius::Float64
    Ball2(center, radius) = radius < 0. ? throw(DomainError()) : new(center, radius)
end

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
function σ(d::Union{Vector{Float64}, SparseVector{Float64,Int64}}, B::Ball2)::Vector{Float64}
    dnorm = norm(d)
    if dnorm > 0
        return B.center .+ d .* (B.radius / dnorm)
    else
        return zeros(length(d))
    end
end

