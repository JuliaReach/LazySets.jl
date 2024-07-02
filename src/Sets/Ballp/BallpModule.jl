module BallpModule

using Reexport

using ..LazySets: AbstractBallp, Ball1, Ball2, BallInf
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Distribution: reseed!

@reexport import ..API: center, isoperationtype, rand, reflect, project, scale,
                        translate!
@reexport import ..LazySets: ball_norm, radius_ball
@reexport using ..API

export Ballp

"""
    Ballp{N<:AbstractFloat, VN<:AbstractVector{N}} <: AbstractBallp{N}

Type that represents a ball in the p-norm, for ``1 ≤ p ≤ ∞``.

It is defined as the set

```math
\\mathcal{B}_p^n(c, r) = \\{ x ∈ ℝ^n : ‖ x - c ‖_p ≤ r \\},
```
where ``c ∈ ℝ^n`` is its center and ``r ∈ ℝ_+`` its radius.
Here ``‖ ⋅ ‖_p`` for ``1 ≤ p ≤ ∞`` denotes the vector ``p``-norm, defined as
``‖ x ‖_p = \\left( ∑\\limits_{i=1}^n |x_i|^p \\right)^{1/p}`` for any
``x ∈ ℝ^n``.

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
Ballp{Float64, Vector{Float64}}(1.5, [0.0, 0.0, 0.0, 0.0, 0.0], 0.5)

julia> dim(B)
5
```

We evaluate the support vector in direction ``[1,2,…,5]``:

```jldoctest ballp_constructor
julia> σ([1.0, 2, 3, 4, 5], B)
5-element Vector{Float64}:
 0.013516004434607558
 0.05406401773843023
 0.12164403991146802
 0.21625607095372093
 0.33790011086518895
```
"""
struct Ballp{N<:AbstractFloat,VN<:AbstractVector{N}} <: AbstractBallp{N}
    p::N
    center::VN
    radius::N

    # default constructor with domain constraint for radius and p
    function Ballp(p::N, center::VN, radius::N) where {N,VN<:AbstractVector{N}}
        @assert radius >= zero(N) "the radius must not be negative"
        @assert p >= one(N) "p must not be less than 1"
        if p == N(Inf)
            return BallInf(center, radius)
        elseif p == N(2)
            return Ball2(center, radius)
        elseif isone(p)
            return Ball1(center, radius)
        else
            return new{N,VN}(p, center, radius)
        end
    end
end

isoperationtype(::Type{<:Ballp}) = false

"""
    center(B::Ballp)

Return the center of a ball in the p-norm.

### Input

- `B` -- ball in the p-norm

### Output

The center of the ball in the p-norm.
"""
function center(B::Ballp)
    return B.center
end

function radius_ball(B::Ballp)
    return B.radius
end

function ball_norm(B::Ballp)
    return B.p
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
              seed::Union{Int,Nothing}=nothing)
    rng = reseed!(rng, seed)
    p = one(N) + abs(randn(rng, N))
    center = randn(rng, N, dim)
    radius = abs(randn(rng, N))
    return Ballp(p, center, radius)
end

"""
    translate!(B::Ballp, v::AbstractVector)

Translate (i.e., shift) a ball in the p-norm by a given vector, in-place.

### Input

- `B` -- ball in the p-norm
- `v` -- translation vector

### Output

The ball `B` translated by `v`.

### Algorithm

We add the vector to the center of the ball.

### Notes

See also [`translate(::Ballp, ::AbstractVector)`](@ref) for the out-of-place
version.
"""
function translate!(B::Ballp, v::AbstractVector)
    @assert length(v) == dim(B) "cannot translate a $(dim(B))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    c = center(B)
    c .+= v
    return B
end

function project(B::Ballp, block::AbstractVector{Int}; kwargs...)
    return Ballp(B.p, B.center[block], B.radius)
end

"""
    reflect(B::Ballp)

Concrete reflection of a ball in the p-norm `B`, resulting in the reflected set
`-B`.

### Input

- `B` -- ball in the p-norm

### Output

The `Ballp` representing `-B`.

### Algorithm

If ``B`` has center ``c`` and radius ``r``, then ``-B`` has center ``-c`` and
radius ``r``. The norm remains the same.
"""
function reflect(B::Ballp)
    return Ballp(B.p, -center(B), B.radius)
end

function scale(α::Real, B::Ballp)
    return Ballp(B.p, B.center .* α, B.radius * abs(α))
end

end  # module
