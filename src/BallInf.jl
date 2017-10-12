export BallInf

"""
    BallInf <: LazySet

Type that represents a ball in the infinity norm,
``B_∞(c, r) = \\{ x ∈ \\mathbb{R}^n : ‖x - c‖_∞ ≦ r \\}.``

### Fields

- ``c`` -- a real vector, the center
- ``r`` -- the radius ( ``≧ 0`` )

### Examples

We create the two-dimensional unit ball, and compute its support function
along the direction `[1, 1]`:

```julia
    julia> B = BallInf(zeros(2), 0.1)
    LazySets.BallInf([0.0, 0.0], 0.1)

    julia> dim(B)
    2

    julia> ρ([1., 1.], B)
    0.2
```
"""
struct BallInf <: LazySet
    center::Vector{Float64}
    radius::Float64
    BallInf(center, radius) = radius < 0. ? throw(DomainError()) : new(center, radius)
end

# dimension of a ball in the infinity norm
function dim(B::BallInf)::Int64
    return length(B.center)
end

"""
    σ(d, P)

Return the support vector of an infinity-norm ball in a given direction.

INPUT:

- ``d`` -- direction
- ``P`` -- polyhedron in H-representation

ALGORITHM:

This code is a vectorized version of
```julia
[(d[i] >= 0) ? B.center[i] + B.radius : B.center[i] - B.radius for i in 1:length(d)]
```
Notice that we cannot use ``B.center + sign.(d) * B.radius``, since the built-in ``sign``
function is such that `sign(0) = 0`, instead of 1. For this reason, we use the
custom ``unit_step`` function, that allows to do: ``B.center + unit_step.(d) * B.radius``
(the dot operator performs broadcasting, to accept vector-valued entries).
"""
function σ(d::Union{Vector{Float64}, SparseVector{Float64,Int64}}, B::BallInf)::Vector{Float64}
    return B.center .+ unit_step.(d) .* B.radius
end

