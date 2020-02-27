import Base: rand,
             ∈

export Ball2,
       sample,
       volume

"""
    Ball2{N<:AbstractFloat} <: AbstractCentrallySymmetric{N}

Type that represents a ball in the 2-norm.

### Fields

- `center` -- center of the ball as a real vector
- `radius` -- radius of the ball as a real scalar (``≥ 0``)

### Notes

Mathematically, a ball in the 2-norm is defined as the set

```math
\\mathcal{B}_2^n(c, r) = \\{ x ∈ \\mathbb{R}^n : ‖ x - c ‖_2 ≤ r \\},
```
where ``c ∈ \\mathbb{R}^n`` is its center and ``r ∈ \\mathbb{R}_+`` its radius.
Here ``‖ ⋅ ‖_2`` denotes the Euclidean norm (also known as 2-norm), defined as
``‖ x ‖_2 = \\left( \\sum\\limits_{i=1}^n |x_i|^2 \\right)^{1/2}`` for any
``x ∈ \\mathbb{R}^n``.

### Examples

Create a five-dimensional ball `B` in the 2-norm centered at the origin with
radius 0.5:

```jldoctest ball2_label
julia> B = Ball2(zeros(5), 0.5)
Ball2{Float64}([0.0, 0.0, 0.0, 0.0, 0.0], 0.5)
julia> dim(B)
5
```

Evaluate `B`'s support vector in the direction ``[1,2,3,4,5]``:

```jldoctest ball2_label
julia> σ([1.,2.,3.,4.,5.], B)
5-element Array{Float64,1}:
 0.06741998624632421
 0.13483997249264842
 0.20225995873897262
 0.26967994498529685
 0.3370999312316211
```
"""
struct Ball2{N<:AbstractFloat} <: AbstractCentrallySymmetric{N}
    center::Vector{N}
    radius::N

    # default constructor with domain constraint for radius
    function Ball2(center::Vector{N}, radius::N) where {N<:AbstractFloat}
        @assert radius >= zero(N) "radius must not be negative"
        return new{N}(center, radius)
    end
end

isoperationtype(::Type{<:Ball2}) = false
isconvextype(::Type{<:Ball2}) = true


# --- AbstractCentrallySymmetric interface functions ---


"""
    center(B::Ball2{N}) where {N<:AbstractFloat}

Return the center of a ball in the 2-norm.

### Input

- `B` -- ball in the 2-norm

### Output

The center of the ball in the 2-norm.
"""
function center(B::Ball2{N}) where {N<:AbstractFloat}
    return B.center
end


# --- LazySet interface functions ---


"""
    σ(d::AbstractVector{N}, B::Ball2{N}) where {N<:AbstractFloat}

Return the support vector of a 2-norm ball in a given direction.

### Input

- `d` -- direction
- `B` -- ball in the 2-norm

### Output

The support vector in the given direction.
If the direction has norm zero, the origin is returned.

### Notes

Let ``c`` and ``r`` be the center and radius of a ball ``B`` in the 2-norm,
respectively.
For nonzero direction ``d`` we have

```math
σ_B(d) = c + r \\frac{d}{‖d‖_2}.
```

This function requires computing the 2-norm of the input direction, which is
performed in the given precision of the numeric datatype of both the direction
and the set.
Exact inputs are not supported.
"""
function σ(d::AbstractVector{N}, B::Ball2{N}) where {N<:AbstractFloat}
    dnorm = norm(d, 2)
    if dnorm <= zero(N)
        return zeros(N, length(d))
    else
        return @. B.center + d * (B.radius / dnorm)
    end
end

"""
    ∈(x::AbstractVector{N}, B::Ball2{N}) where {N<:AbstractFloat}

Check whether a given point is contained in a ball in the 2-norm.

### Input

- `x` -- point/vector
- `B` -- ball in the 2-norm

### Output

`true` iff ``x ∈ B``.

### Notes

This implementation is worst-case optimized, i.e., it is optimistic and first
computes (see below) the whole sum before comparing to the radius.
In applications where the point is typically far away from the ball, a fail-fast
implementation with interleaved comparisons could be more efficient.

### Algorithm

Let ``B`` be an ``n``-dimensional ball in the 2-norm with radius ``r`` and let
``c_i`` and ``x_i`` be the ball's center and the vector ``x`` in dimension
``i``, respectively.
Then ``x ∈ B`` iff ``\\left( ∑_{i=1}^n |c_i - x_i|^2 \\right)^{1/2} ≤ r``.

### Examples

```jldoctest
julia> B = Ball2([1., 1.], sqrt(0.5))
Ball2{Float64}([1.0, 1.0], 0.7071067811865476)
julia> [.5, 1.6] ∈ B
false
julia> [.5, 1.5] ∈ B
true
```
"""
function ∈(x::AbstractVector{N}, B::Ball2{N}) where {N<:AbstractFloat}
    @assert length(x) == dim(B)
    sum = zero(N)
    @inbounds for i in eachindex(x)
        sum += (B.center[i] - x[i])^2
    end
    return _leq(sqrt(sum), B.radius)
end

"""
    rand(::Type{Ball2}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

Create a random ball in the 2-norm.

### Input

- `Ball2` -- type for dispatch
- `N`     -- (optional, default: `Float64`) numeric type
- `dim`   -- (optional, default: 2) dimension
- `rng`   -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`  -- (optional, default: `nothing`) seed for reseeding

### Output

A random ball in the 2-norm.

### Algorithm

All numbers are normally distributed with mean 0 and standard deviation 1.
Additionally, the radius is nonnegative.
"""
function rand(::Type{Ball2};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int, Nothing}=nothing)
    rng = reseed(rng, seed)
    center = randn(rng, N, dim)
    radius = abs(randn(rng, N))
    return Ball2(center, radius)
end

"""
    translate(B::Ball2{N}, v::AbstractVector{N}) where {N<:AbstractFloat}

Translate (i.e., shift) a ball in the 2-norm by a given vector.

### Input

- `B` -- ball in the 2-norm
- `v` -- translation vector

### Output

A translated ball in the 2-norm.

### Algorithm

We add the vector to the center of the ball.
"""
function translate(B::Ball2{N}, v::AbstractVector{N}) where {N<:AbstractFloat}
    @assert length(v) == dim(B) "cannot translate a $(dim(B))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    return Ball2(center(B) + v, B.radius)
end

"""
    sample(B::Ball2{N}, nsamples::Int=1;
           [rng]::AbstractRNG=GLOBAL_RNG,
           [seed]::Union{Int, Nothing}=nothing) where {N<:AbstractFloat}

Return samples from a uniform distribution on the given ball in the 2-norm.

### Input

- `B`        -- ball in the 2-norm
- `nsamples` -- (optional, default: `1`) number of samples
- `rng`      -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`     -- (optional, default: `nothing`) seed for reseeding

### Output

A linear array of `nsamples` elements drawn from a uniform distribution in `B`.

### Algorithm

Random sampling with uniform distribution in `B` is computed using Muller's method
of normalized Gaussians. This function requires the package `Distributions`.
See `_sample_unit_nball_muller!` for implementation details.
"""
function sample(B::Ball2{N}, nsamples::Int=1;
                rng::AbstractRNG=GLOBAL_RNG,
                seed::Union{Int, Nothing}=nothing) where {N<:AbstractFloat}
    require(:Distributions; fun_name="sample")
    n = dim(B)
    D = Vector{Vector{N}}(undef, nsamples) # preallocate output
    _sample_unit_nball_muller!(D, n, nsamples, rng=rng, seed=seed)

    # customize for the given ball
    r, c = B.radius, B.center
    @inbounds for i in 1:nsamples
        axpby!(one(N), c, r, D[i])
    end
    return D
end


# --- Ball2 functions ---


"""
    chebyshev_center(B::Ball2{N}; compute_radius::Bool=false
                    ) where {N<:AbstractFloat}

Compute the [Chebyshev center](https://en.wikipedia.org/wiki/Chebyshev_center)
of a ball in the 2-norm.

### Input

- `B`              -- ball in the 2-norm
- `compute_radius` -- (optional; default: `false`) option to additionally return
                      the radius of the largest ball enclosed by `B` around the
                      Chebyshev center

### Output

If `compute_radius` is `false`, the result is the Chebyshev center of `B`.
If `compute_radius` is `true`, the result is the pair `(c, r)` where `c` is the
Chebyshev center of `B` and `r` is the radius of the largest ball with center
`c` enclosed by `B`.

### Notes

The Chebyshev center of a ball in the 2-norm is just the center of the ball.
"""
function chebyshev_center(B::Ball2{N}; compute_radius::Bool=false
                         ) where {N<:AbstractFloat}
    if compute_radius
        return center(B), B.radius
    end
    return center(B)
end

"""
    volume(B::Ball2{N}) where {N<:AbstractFloat}

Return the volume of a ball in the 2-norm.

### Input

- `B` -- ball in the 2-norm

### Output

The volume of ``B``.

### Algorithm

This function implements the well-known formula for the volume of an n-dimensional
ball using factorials. For details see the wikipedia article
[Volume of an n-ball](https://en.wikipedia.org/wiki/Volume_of_an_n-ball).
"""
function volume(B::Ball2{N}) where {N<:AbstractFloat}
    n = dim(B)
    k = div(n, 2)
    R = radius(B)
    if iseven(n)
        vol = Base.pi^k * R^n / factorial(k)
    else
        vol = 2 * factorial(k) * (4*Base.pi)^k * R^n / factorial(n)
    end
    return vol
end
