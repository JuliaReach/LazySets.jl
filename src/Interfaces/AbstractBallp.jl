export AbstractBallp

"""
    AbstractBallp{N} <: AbstractCentrallySymmetric{N}

Abstract type for p-norm balls.

### Notes

See [`Ballp`](@ref) for a standard implementation of this interface.

Every concrete `AbstractBallp` must define the following methods:

- `ball_norm(::AbstractBallp)` -- return the characteristic norm
- `radius_ball(::AbstractBallp)` -- return the ball radius

The subtypes of `AbstractBallp`:

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(AbstractBallp)
2-element Vector{Any}:
 Ball2
 Ballp
```

There are two further set types implementing the `AbstractBallp` interface, but
they also implement other interfaces and hence cannot be subtypes: `Ball1` and
`BallInf`.
"""
abstract type AbstractBallp{N} <: AbstractCentrallySymmetric{N} end

"""
    radius_ball(B::AbstractBallp)

Compute the radius of a p-norm ball.

### Input

- `B` -- p-norm ball

### Output

A number representing the radius.
"""
function radius_ball(::AbstractBallp) end

"""
    ball_norm(B::AbstractBallp)

Determine the norm (p) of a p-norm ball.

### Input

- `B` -- p-norm ball

### Output

A number representing the norm.
"""
function ball_norm(::AbstractBallp) end

function low(B::AbstractBallp)
    return _low_AbstractBallp(B)
end

function _low_AbstractBallp(B::LazySet)
    return center(B) .- radius_ball(B)
end

function low(B::AbstractBallp, i::Int)
    return _low_AbstractBallp(B, i)
end

function _low_AbstractBallp(B::LazySet, i::Int)
    return center(B, i) - radius_ball(B)
end

function high(B::AbstractBallp)
    return _high_AbstractBallp(B)
end

function _high_AbstractBallp(B::LazySet)
    return center(B) .+ radius_ball(B)
end

function high(B::AbstractBallp, i::Int)
    return _high_AbstractBallp(B, i)
end

function _high_AbstractBallp(B::LazySet, i::Int)
    return center(B, i) + radius_ball(B)
end

"""
# Extended help

    σ(d::AbstractVector, B::AbstractBallp)

### Algorithm

The support vector of the unit ball in the ``p``-norm along direction ``d`` is:

```math
σ(d, \\mathcal{B}_p^n(0, 1)) = \\dfrac{\\tilde{v}}{‖\\tilde{v}‖_q},
```
where ``\\tilde{v}_i = \\frac{|d_i|^q}{d_i}`` if ``d_i ≠ 0`` and
``\\tilde{v}_i = 0`` otherwise, for all ``i=1,…,n``, and ``q`` is the conjugate
number of ``p``.
By the affine transformation ``x = r\\tilde{x} + c``, one obtains that
the support vector of ``\\mathcal{B}_p^n(c, r)`` is

```math
σ(d, \\mathcal{B}_p^n(c, r)) = \\dfrac{v}{‖v‖_q},
```
where ``v_i = c_i + r\\frac{|d_i|^q}{d_i}`` if ``d_i ≠ 0`` and ``v_i = 0``
otherwise, for all ``i = 1, …, n``.

If the direction has norm zero, the center of the ball is returned.
"""
@validate function σ(d::AbstractVector, B::AbstractBallp)
    p = ball_norm(B)
    q = p / (p - 1)
    v = similar(d)
    N = promote_type(eltype(d), eltype(B))
    @inbounds for (i, di) in enumerate(d)
        v[i] = di == zero(N) ? di : abs.(di) .^ q / di
    end
    vnorm = norm(v, p)
    if isapproxzero(vnorm)
        svec = copy(center(B))
    else
        svec = center(B) .+ radius_ball(B) .* (v ./ vnorm)
    end
    return svec
end

"""
# Extended help

    ρ(d::AbstractVector, B::AbstractBallp)

### Algorithm

Let ``c`` and ``r`` be the center and radius of the ball ``B`` in the p-norm,
respectively, and let ``q = \\frac{p}{p-1}``. Then:

```math
ρ(d, B) = ⟨d, c⟩ + r ‖d‖_q.
```
"""
@validate function ρ(d::AbstractVector, B::AbstractBallp)
    p = ball_norm(B)
    q = p / (p - 1)
    return dot(d, center(B)) + radius_ball(B) * norm(d, q)
end

"""
# Extended help

    ∈(x::AbstractVector, B::AbstractBallp)

### Notes

This implementation is worst-case optimized, i.e., it is optimistic and first
computes (see below) the whole sum before comparing to the radius.
In applications where the point is typically far away from the ball, a fail-fast
implementation with interleaved comparisons could be more efficient.

### Algorithm

Let ``B`` be an ``n``-dimensional ball in the p-norm with radius ``r`` and let
``c_i`` and ``x_i`` be the ball's center and the vector ``x`` in dimension
``i``, respectively.
Then ``x ∈ B`` iff ``\\left( ∑_{i=1}^n |c_i - x_i|^p \\right)^{1/p} ≤ r``.

### Examples

```jldoctest
julia> B = Ballp(1.5, [1.0, 1.0], 1.)
Ballp{Float64, Vector{Float64}}(1.5, [1.0, 1.0], 1.0)

julia> [0.5, -0.5] ∈ B
false

julia> [0.5, 1.5] ∈ B
true
```
"""
@validate function ∈(x::AbstractVector, B::AbstractBallp)
    N = promote_type(eltype(x), eltype(B))
    p = ball_norm(B)
    sum = zero(N)
    @inbounds for i in eachindex(x)
        sum += abs(center(B, i) - x[i])^p
    end
    return _leq(sum^(one(N) / p), radius_ball(B))
end
