"""
    Ball2 <: LazySet

Type that represents a ball in the 2-norm.

FIELDS:

- ``c`` -- a real vector, the center
- ``r`` -- the radius (>= 0)

EXAMPLES:

A ten-dimensional ball in the 2-norm centered at the origin, and of radius 0.5:

    julia> B = Ball2(zeros(10), 0.5)
    Ball2([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 0.5)
    julia> dim(B)
    10
    julia> σ(collect(range(0, 0.1, 10)), B)
    10-element Array{Float64,1}:
     0.0
     0.0296174
     0.0592349
     0.0888523
     0.11847
     0.148087
     0.177705
     0.207322
     0.23694
     0.266557
"""
struct Ball2 <: LazySet
    center::Vector{Float64}
    radius::Float64
    Ball2(center, radius) = radius < 0. ? throw(DomainError()) : new(center, radius)
end

"""
    dim(H)

Return the dimension of a Ball2.

INPUT:

- ``B`` -- a ball in the 2-norm

OUTPUT:

The ambient dimension of the ball.
"""
function dim(B::Ball2)::Int64
    length(B.center)
end

r"""
    σ(d, B)

Return the support vector of a Ball2 in a given direction.

INPUT:

- ``d`` -- a direction
- ``B`` -- a ball in the 2-norm

OUTPUT:

The support vector in the given direction.

NOTES:

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

export Ball2
