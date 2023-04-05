import Base.rand

export BallInf,
       volume

const BALLINF_THRESHOLD_ρ = 30  # threshold value in `ρ`
const BALLINF_THRESHOLD_VOLUME = 50  # threshold value in `volume`

"""
    BallInf{N, VN<:AbstractVector{N}} <: AbstractHyperrectangle{N}

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

Construct the two-dimensional unit ball and compute its support function along
the positive ``x=y`` direction:

```jldoctest
julia> B = BallInf(zeros(2), 1.0)
BallInf{Float64, Vector{Float64}}([0.0, 0.0], 1.0)

julia> dim(B)
2

julia> ρ([1.0, 1.0], B)
2.0
```
"""
struct BallInf{N, VN<:AbstractVector{N}} <: AbstractHyperrectangle{N}
    center::VN
    radius::N

    # default constructor with domain constraint for radius
    function BallInf(center::VN, radius::N) where {N, VN<:AbstractVector{N}}
        @assert radius >= zero(N) "the radius must not be negative"
        return new{N, VN}(center, radius)
    end
end

isoperationtype(::Type{<:BallInf}) = false

"""
    radius_hyperrectangle(B::BallInf, i::Int)

Return the box radius of a ball in the infinity norm in a given dimension.

### Input

- `B` -- ball in the infinity norm
- `i` -- dimension of interest

### Output

The box radius of the ball in the infinity norm in the given dimension.
"""
function radius_hyperrectangle(B::BallInf, i::Int)
    @assert 1 <= i <= dim(B) "cannot compute the radius of a " *
        "$(dim(B))-dimensional set in dimension $i"
    return B.radius
end

"""
    radius_hyperrectangle(B::BallInf)

Return the box radius of a ball in the infinity norm.

### Input

- `B` -- ball in the infinity norm

### Output

The box radius of the ball in the infinity norm, which is the same in every
dimension.
"""
function radius_hyperrectangle(B::BallInf)
    return fill(B.radius, dim(B))
end

"""
    isflat(B::BallInf)

Determine whether a ball in the infinity norm is flat, i.e., whether its radius
is zero.

### Input

- `B` -- ball in the infinity norm

### Output

`true` iff the ball is flat.

### Notes

For robustness with respect to floating-point inputs, this function relies on
the result of `isapproxzero` applied to the radius of the ball.
Hence, this function depends on the absolute zero tolerance `ABSZTOL`.
"""
function isflat(B::BallInf)
    return isapproxzero(B.radius)
end

function load_genmat_ballinf_static()
return quote

function genmat(B::BallInf{N, SVector{L, N}}) where {L, N}
    if isflat(B)
        return SMatrix{L, 0, N, 0}()
    else
        gens = zeros(MMatrix{L, L})
        @inbounds for i in 1:L
            gens[i, i] = B.radius
        end
        return SMatrix(gens)
    end
end

end end  # quote / load_genmat_ballinf_static()

"""
    center(B::BallInf)

Return the center of a ball in the infinity norm.

### Input

- `B` -- ball in the infinity norm

### Output

The center of the ball in the infinity norm.
"""
function center(B::BallInf)
    return B.center
end

"""
    σ(d::AbstractVector, B::BallInf)

Return the support vector of a ball in the infinity norm in the given direction.

### Input

- `d` -- direction
- `B` -- ball in the infinity norm

### Output

The support vector in the given direction.
If the direction has norm zero, the center of the ball is returned.
"""
function σ(d::AbstractVector, B::BallInf)
    @assert length(d) == dim(B) "a $(length(d))-dimensional vector is " *
                                "incompatible with a $(dim(B))-dimensional set"
    return center(B) .+ sign.(d) .* B.radius
end

# special case for SingleEntryVector
function σ(d::SingleEntryVector, B::BallInf)
    return _σ_sev_hyperrectangle(d, B)
end

"""
    ρ(d::AbstractVector, B::BallInf)

Evaluate the support function of a ball in the infinity norm in the given
direction.

### Input

- `d` -- direction
- `B` -- ball in the infinity norm

### Output

Evaluation of the support function in the given direction.

### Algorithm

Let ``B`` be a ball in the infinity norm with center ``c`` and radius ``r`` and
let ``d`` be the direction of interest.
For balls with dimensions less than 30 we use the implementation for
`AbstractHyperrectangle`, taylored to a `BallInf`, which computes

```math
    ∑_{i=1}^n d_i (c_i + \\textrm{sgn}(d_i) · r)
```

where ``\\textrm{sgn}(α) = 1`` if ``α ≥ 0`` and ``\\textrm{sgn}(α) = -1`` if ``α < 0``.

For balls of higher dimension we instead exploit that for a support vector
``v = σ(d, B) = c + \\textrm{sgn}(d) · (r, …, r)ᵀ`` we have

```math
    ρ(d, B) = ⟨d, v⟩ = ⟨d, c⟩ + ⟨d, \\textrm{sgn}(d) · (r, …, r)ᵀ⟩ = ⟨d, c⟩ + r · ∑_{i=1}^n |d_i|
```

where ``⟨·, ·⟩`` denotes the dot product.
"""
function ρ(d::AbstractVector, B::BallInf)
    @assert length(d) == dim(B) "a $(length(d))-dimensional vector is " *
                                "incompatible with a $(dim(B))-dimensional set"
    c = center(B)
    r = B.radius
    if length(d) > BALLINF_THRESHOLD_ρ
        # more efficient for higher dimensions
        return dot(d, c) + r * sum(abs, d)
    end
    N = promote_type(eltype(d), eltype(B))
    res = zero(N)
    @inbounds for (i, di) in enumerate(d)
        if di < zero(N)
            res += di * (c[i] - r)
        elseif di > zero(N)
            res += di * (c[i] + r)
        end
    end
    return res
end

# special case for SingleEntryVector
function ρ(d::SingleEntryVector, B::BallInf)
    return _ρ_sev_hyperrectangle(d, B)
end

"""
    radius(B::BallInf, [p]::Real=Inf)

Return the radius of a ball in the infinity norm.

### Input

- `B` -- ball in the infinity norm
- `p` -- (optional, default: `Inf`) norm

### Output

A real number representing the radius.

### Notes

The result is defined as the radius of the enclosing ball of the given
``p``-norm of minimal volume with the same center.
"""
function radius(B::BallInf, p::Real=Inf)
    return (p == Inf) ? B.radius : norm(fill(B.radius, dim(B)), p)
end

"""
    rand(::Type{BallInf}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

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
              seed::Union{Int, Nothing}=nothing)
    rng = reseed(rng, seed)
    center = randn(rng, N, dim)
    radius = abs(randn(rng, N))
    return BallInf(center, radius)
end

"""
    translate(B::BallInf, v::AbstractVector)

Translate (i.e., shift) a ball in the infinity norm by a given vector.

### Input

- `B` -- ball in the infinity norm
- `v` -- translation vector

### Output

A translated ball in the infinity norm.

### Notes

See also [`translate!(::BallInf, ::AbstractVector)`](@ref) for the in-place
version.
"""
function translate(B::BallInf, v::AbstractVector)
    return translate!(copy(B), v)
end

"""
    translate!(B::BallInf, v::AbstractVector)

Translate (i.e., shift) a ball in the infinity norm by a given vector, in-place.

### Input

- `B` -- ball in the infinity norm
- `v` -- translation vector

### Output

The ball `B` translated by `v`.

### Algorithm

We add the vector to the center of the ball.

### Notes

See also [`translate(::BallInf, ::AbstractVector)`](@ref) for the out-of-place
version.
"""
function translate!(B::BallInf, v::AbstractVector)
    @assert length(v) == dim(B) "cannot translate a $(dim(B))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    c = B.center
    c .+= v
    return B
end

# compute a^n in a loop
@inline function _pow_loop(a::N, n::Int) where {N}
    vol = one(N)
    diam = 2 * a
    for i in 1:n
        vol *= diam
    end
    return vol
end

# compute a^n using exp
@inline function _pow_exp(a, n::Int)
    return exp(n * log(2a))
end

"""
    volume(B::BallInf)

Return the volume of a ball in the infinity norm.

### Input

- `B` -- ball in the infinity norm

### Output

The volume of ``B``.

### Algorithm

We compute the volume by iterative multiplication of the radius.

For floating-point inputs we use this implementation for balls of dimension less
than 50. For balls of higher dimension we instead compute ``exp(n * log(2r))``,
where ``r`` is the radius of the ball.
"""
function volume(B::BallInf)
    return _pow_loop(B.radius, dim(B))
end

# method for floating-point input
function volume(B::BallInf{N}) where {N<:AbstractFloat}
    n = dim(B)
    if n < BALLINF_THRESHOLD_VOLUME
        vol = _pow_loop(B.radius, n)
    else
        vol = _pow_exp(B.radius, n)
    end
    return vol
end

"""
    ngens(B::BallInf)

Return the number of generators of a ball in the infinity norm.

### Input

- `B` -- ball in the infinity norm

### Output

The number of generators.

### Algorithm

A ball in the infinity norm has either one generator for each dimension, or zero
generators if it is a degenerated ball of radius zero.
"""
function ngens(B::BallInf)
    return iszero(B.radius) ? 0 : dim(B)
end

function project(B::BallInf, block::AbstractVector{Int}; kwargs...)
    return BallInf(B.center[block], B.radius)
end

"""
    reflect(B::BallInf)

Concrete reflection of a ball in the infinity norm `B`, resulting in the
reflected set `-B`.

### Input

- `B` -- ball in the infinity norm

### Output

The `BallInf` representing `-B`.

### Algorithm

If ``B`` has center ``c`` and radius ``r``, then ``-B`` has center ``-c`` and
radius ``r``.
"""
function reflect(B::BallInf)
    return BallInf(-center(B), B.radius)
end
