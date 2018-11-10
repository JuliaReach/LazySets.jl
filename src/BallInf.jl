import Base.rand

export BallInf

"""
    BallInf{N<:Real} <: AbstractHyperrectangle{N}

Type that represents a ball in the infinity norm.

### Fields

- `center` -- center of the ball as a real vector
- `radius` -- radius of the ball as a real scalar (``≥ 0``)

### Notes

Mathematically, a ball in the infinity norm is defined as the set

```math
\\mathcal{B}_∞^n(c, r) = \\{ x ∈ \\mathbb{R}^n : ‖ x - c ‖_∞ ≤ r \\},
```
where ``c ∈ \\mathbb{R}^n`` is its center and ``r ∈ \\mathbb{R}_+`` its radius.
Here ``‖ ⋅ ‖_∞`` denotes the infinity norm, defined as
``‖ x ‖_∞ = \\max\\limits_{i=1,…,n} \\vert x_i \\vert`` for any
``x ∈ \\mathbb{R}^n``.

### Examples

Create the two-dimensional unit ball and compute its support function along the
positive ``x=y`` direction:

```jldoctest
julia> B = BallInf(zeros(2), 1.0)
BallInf{Float64}([0.0, 0.0], 1.0)
julia> dim(B)
2
julia> ρ([1., 1.], B)
2.0
```
"""
struct BallInf{N<:Real} <: AbstractHyperrectangle{N}
    center::Vector{N}
    radius::N

    # default constructor with domain constraint for radius
    function BallInf{N}(center::Vector{N}, radius::N) where {N<:Real}
        @assert radius >= zero(N) "radius must not be negative"
        return new{N}(center, radius)
    end
end

# convenience constructor without type parameter
BallInf(center::Vector{N}, radius::N) where {N<:Real} =
    BallInf{N}(center, radius)


# --- AbstractHyperrectangle interface functions ---


"""
    radius_hyperrectangle(B::BallInf{N}, i::Int)::N where {N<:Real}

Return the box radius of a infinity norm ball in a given dimension.

### Input

- `B` -- infinity norm ball

### Output

The box radius of the ball in the infinity norm in the given dimension.
"""
function radius_hyperrectangle(B::BallInf{N}, i::Int)::N where {N<:Real}
    return B.radius
end

"""
    radius_hyperrectangle(B::BallInf{N})::Vector{N} where {N<:Real}

Return the box radius of a infinity norm ball, which is the same in every
dimension.

### Input

- `B` -- infinity norm ball

### Output

The box radius of the ball in the infinity norm.
"""
function radius_hyperrectangle(B::BallInf{N})::Vector{N} where {N<:Real}
    return fill(B.radius, dim(B))
end


# --- AbstractCentrallySymmetric interface functions ---


"""
    center(B::BallInf{N})::Vector{N} where {N<:Real}

Return the center of a ball in the infinity norm.

### Input

- `B` -- ball in the infinity norm

### Output

The center of the ball in the infinity norm.
"""
function center(B::BallInf{N})::Vector{N} where {N<:Real}
    return B.center
end


# --- LazySet interface functions ---


"""
    radius(B::BallInf, [p]::Real=Inf)::Real

Return the radius of a ball in the infinity norm.

### Input

- `B` -- ball in the infinity norm
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the radius.

### Notes

The radius is defined as the radius of the enclosing ball of the given
``p``-norm of minimal volume with the same center.
"""
function radius(B::BallInf, p::Real=Inf)::Real
    return (p == Inf) ? B.radius : norm(fill(B.radius, dim(B)), p)
end

"""
    rand(::Type{BallInf}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing
        )::BallInf{N}

Create a random ball in the infinity norm.

### Input

- `BallInf` -- type for dispatch
- `N`       -- (optional, default: `Float64`) numeric type
- `dim`     -- (optional, default: 2) dimension
- `rng`     -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`    -- (optional, default: `nothing`) seed for reseeding

### Output

A random ball in the infinity norm.

### Algorithm

All numbers are normally distributed with mean 0 and standard deviation 1.
Additionally, the radius is nonnegative.
"""
function rand(::Type{BallInf};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int, Nothing}=nothing
             )::BallInf{N}
    rng = reseed(rng, seed)
    center = randn(rng, N, dim)
    radius = abs(randn(rng, N))
    return BallInf(center, radius)
end
