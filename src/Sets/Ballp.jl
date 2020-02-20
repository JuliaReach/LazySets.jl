import Base: rand,
             ∈

export Ballp

"""
    Ballp{N<:AbstractFloat} <: AbstractCentrallySymmetric{N}

Type that represents a ball in the p-norm, for ``1 ≤ p ≤ ∞``.

It is defined as the set

```math
\\mathcal{B}_p^n(c, r) = \\{ x ∈ \\mathbb{R}^n : ‖ x - c ‖_p ≤ r \\},
```
where ``c ∈ \\mathbb{R}^n`` is its center and ``r ∈ \\mathbb{R}_+`` its radius.
Here ``‖ ⋅ ‖_p`` for ``1 ≤ p ≤ ∞`` denotes the vector ``p``-norm, defined as
``‖ x ‖_p = \\left( \\sum\\limits_{i=1}^n |x_i|^p \\right)^{1/p}`` for any
``x ∈ \\mathbb{R}^n``.

### Fields

- `p`      -- norm as a real scalar
- `center` -- center of the ball as a real vector
- `radius` -- radius of the ball as a scalar (``≥ 0``)

## Notes

The special cases ``p=1``, ``p=2`` and ``p=∞`` fall back to the specialized
types `Ball1`, `Ball2` and `BallInf`, respectively.

### Examples

A five-dimensional ball in the ``p=3/2`` norm centered at the origin of radius
0.5:

```jldoctest ballp_constructor
julia> B = Ballp(3/2, zeros(5), 0.5)
Ballp{Float64}(1.5, [0.0, 0.0, 0.0, 0.0, 0.0], 0.5)
julia> dim(B)
5
```

We evaluate the support vector in direction ``[1,2,…,5]``:

```jldoctest ballp_constructor
julia> σ([1., 2, 3, 4, 5], B)
5-element Array{Float64,1}:
 0.013516004434607558
 0.05406401773843023
 0.12164403991146802
 0.21625607095372093
 0.33790011086518895
```
"""
struct Ballp{N<:AbstractFloat} <: AbstractCentrallySymmetric{N}
    p::N
    center::Vector{N}
    radius::N

    # default constructor with domain constraint for radius and p
    function Ballp{N}(p::N, center::Vector{N}, radius::N
                     ) where {N<:Real}
        @assert radius >= zero(N) "radius must not be negative"
        @assert p >= 1 "p must not be less than 1"
        if p == Inf
            return BallInf(center, radius)
        elseif p == 2
            return Ball2(center, radius)
        elseif p == 1
            return Ball1(center, radius)
        else
            return new{N}(p, center, radius)
        end
    end
end

isoperationtype(::Type{<:Ballp}) = false
isconvextype(::Type{<:Ballp}) = true

# convenience constructor without type parameter
function Ballp(p::N, center::Vector{N}, radius::N) where {N<:Real}
    if p == Inf
        return BallInf(center, radius)
    elseif p == 2
        return Ball2(center, radius)
    elseif p == 1
        return Ball1(center, radius)
    else
        return Ballp{N}(p, center, radius)
    end
end


# --- AbstractCentrallySymmetric interface functions ---


"""
    center(B::Ballp{N}) where {N<:AbstractFloat}

Return the center of a ball in the p-norm.

### Input

- `B` -- ball in the p-norm

### Output

The center of the ball in the p-norm.
"""
function center(B::Ballp{N}) where {N<:AbstractFloat}
    return B.center
end


# --- LazySet interface functions ---


"""
    σ(d::AbstractVector{N}, B::Ballp{N}) where {N<:AbstractFloat}

Return the support vector of a `Ballp` in a given direction.

### Input

- `d` -- direction
- `B` -- ball in the p-norm

### Output

The support vector in the given direction.
If the direction has norm zero, the center of the ball is returned.

### Algorithm

The support vector of the unit ball in the ``p``-norm along direction ``d`` is:

```math
σ_{\\mathcal{B}_p^n(0, 1)}(d) = \\dfrac{\\tilde{v}}{‖\\tilde{v}‖_q},
```
where ``\\tilde{v}_i = \\frac{|d_i|^q}{d_i}`` if ``d_i ≠ 0`` and
``tilde{v}_i = 0`` otherwise, for all ``i=1,…,n``, and ``q`` is the conjugate
number of ``p``.
By the affine transformation ``x = r\\tilde{x} + c``, one obtains that
the support vector of ``\\mathcal{B}_p^n(c, r)`` is

```math
σ_{\\mathcal{B}_p^n(c, r)}(d) = \\dfrac{v}{‖v‖_q},
```
where ``v_i = c_i + r\\frac{|d_i|^q}{d_i}`` if ``d_i ≠ 0`` and ``v_i = 0``
otherwise, for all ``i = 1, …, n``.
"""
function σ(d::AbstractVector{N}, B::Ballp{N}) where {N<:AbstractFloat}
    p = B.p
    q = p/(p-1)
    v = similar(d)
    @inbounds for (i, di) in enumerate(d)
        v[i] = di == zero(N) ? di : abs.(di).^q / di
    end
    vnorm = norm(v, p)
    svec = vnorm != zero(N) ? @.(B.center + B.radius * (v/vnorm)) : B.center
    return svec
end

"""
    ∈(x::AbstractVector{N}, B::Ballp{N}) where {N<:AbstractFloat}

Check whether a given point is contained in a ball in the p-norm.

### Input

- `x` -- point/vector
- `B` -- ball in the p-norm

### Output

`true` iff ``x ∈ B``.

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
julia> B = Ballp(1.5, [1., 1.], 1.)
Ballp{Float64}(1.5, [1.0, 1.0], 1.0)
julia> [.5, -.5] ∈ B
false
julia> [.5, 1.5] ∈ B
true
```
"""
function ∈(x::AbstractVector{N}, B::Ballp{N}) where {N<:AbstractFloat}
    @assert length(x) == dim(B)
    sum = zero(N)
    for i in eachindex(x)
        sum += abs(B.center[i] - x[i])^B.p
    end
    return sum^(one(N)/B.p) <= B.radius
end

"""
    rand(::Type{Ballp}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

Create a random ball in the p-norm.

### Input

- `Ballp` -- type for dispatch
- `N`     -- (optional, default: `Float64`) numeric type
- `dim`   -- (optional, default: 2) dimension
- `rng`   -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`  -- (optional, default: `nothing`) seed for reseeding

### Output

A random ball in the p-norm.

### Algorithm

The center and radius are normally distributed with mean 0 and standard
deviation 1.
Additionally, the radius is nonnegative.
The p-norm is a normally distributed number ≥ 1 with mean 1 and standard
deviation 1.
"""
function rand(::Type{Ballp};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int, Nothing}=nothing)
    rng = reseed(rng, seed)
    p = one(N) + abs(randn(rng, N))
    center = randn(rng, N, dim)
    radius = abs(randn(rng, N))
    return Ballp(p, center, radius)
end

"""
    translate(B::Ballp{N}, v::AbstractVector{N}) where {N<:AbstractFloat}

Translate (i.e., shift) a ball in the p-norm by a given vector.

### Input

- `B` -- ball in the p-norm
- `v` -- translation vector

### Output

A translated ball in the p- norm.

### Algorithm

We add the vector to the center of the ball.
"""
function translate(B::Ballp{N}, v::AbstractVector{N}) where {N<:AbstractFloat}
    return translate!(copy(B), v)
end

function translate!(B::Ballp{N}, v::AbstractVector{N}) where {N<:AbstractFloat}
    @assert length(v) == dim(B) "cannot translate a $(dim(B))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    c = center(B)
    c .+= v
    return B
end
