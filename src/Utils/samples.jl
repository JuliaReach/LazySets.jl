# ========================
# Sampling from a LazySet
# ========================

"""
    Sampler

Abstract type for defining new sample methods.

### Notes

All subtypes should implement a `_sample!` method.
"""
abstract type Sampler end

"""
    sample(X::LazySet{N}, num_samples::Int;
           [sampler]=_default_sampler(X),
           [rng]::AbstractRNG=GLOBAL_RNG,
           [seed]::Union{Int, Nothing}=nothing,
           [VN]=Vector{N}) where {N}

Sampling of an arbitrary bounded set `X`.

### Input

- `X`           -- (bounded) set to be sampled
- `num_samples` -- number of random samples
- `sampler`     -- (optional, default: `_default_sampler(X)`) the sampler used;
                   falls back to `RejectionSampler` or `UniformSampler` depending
                   on the type of `X` 
- `rng`         -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`        -- (optional, default: `nothing`) seed for reseeding
- `VN`          -- (optional, default: `Vector{N}`) vector type of the sampled points

### Output

A vector of `num_samples` vectors.
If `num_samples` is not passed, the result is just one sample (not wrapped in a
vector).

### Algorithm

See the documentation of the respective `Sampler`.
"""
function sample(X::LazySet{N}, num_samples::Int;
                sampler=_default_sampler(X),
                rng::AbstractRNG=GLOBAL_RNG,
                seed::Union{Int, Nothing}=nothing,
                VN=Vector{N}) where {N}
    @assert isbounded(X) "this function requires that the set `X` is bounded"

    D = Vector{VN}(undef, num_samples) # preallocate output
    _sample!(D, sampler(X); rng=rng, seed=seed)
    return D
end

# without argument, returns a single element (instead of a singleton)
function sample(X::LazySet{N}; kwargs...) where {N}
    return sample(X, 1; kwargs...)[1]
end

# fallback implementation
function _sample!(D::Vector{VN},
                  sampler::Sampler;
                  rng::AbstractRNG=GLOBAL_RNG,
                  seed::Union{Int, Nothing}=nothing) where {N, VN<:AbstractVector{N}}
    error("the method `_sample!` is not implemented for samplers of type " *
          "$(typeof(sampler))")
end

# =====================
# Rejection Sampling
# =====================

# represents a uniform distribution over the interval [a, b]
struct DefaultUniform{N}
    a::N
    b::N
end

function Base.rand(rng::AbstractRNG, U::DefaultUniform)
    (U.b - U.a) * rand(rng) + U.a
end

"""
    RejectionSampler{S<:LazySet, D} <: Sampler

Type used for rejection sampling of an arbitrary `LazySet` `X`.

### Fields

- `X`          -- (bounded) set to be sampled
- `box_approx` -- Distribution from which the sample is drawn

### Algorithm

Draw a sample ``x`` from a uniform distribution of a box-overapproximation of the
original set ``X`` in all ``n`` dimensions. The function rejects a drawn sample ``x``
and redraws as long as the sample is not contained in the original set ``X``,
i.e., ``x ∉ X``.
"""
struct RejectionSampler{S<:LazySet, D} <: Sampler
    X::S
    box_approx::Vector{D}
end

function RejectionSampler(X, distribution=DefaultUniform)
    B = box_approximation(X)
    box_approx = [distribution(low(B, i), high(B, i)) for i in 1:dim(B)]
    return RejectionSampler(X, box_approx)
end

set(sampler::RejectionSampler) = sampler.X

"""
    UniformSampler{S<:LazySet, D} <: Sampler

Type used for uniform sampling of an arbitrary `LazySet` `X`.

### Fields

- `X` -- set to be sampled

### Algorithm

Draw a sample ``x`` from a uniform distribution of the original set ``X``.
"""
struct UniformSampler{S<:LazySet, D} <: Sampler
    X::S
end

UniformSampler(X::LazySet) = UniformSampler{typeof(X), DefaultUniform{eltype(X)}}(X)

set(sampler::UniformSampler) = sampler.X

"""
    PolytopeSampler{S<:LazySet, D} <: Sampler

Type used for sampling of a convex polytope `X`.

### Fields

- `X` -- convex polytope to be sampled

### Algorithm

Choose a random convex combination of the vertices of `X`.

Let ``V = \\{v_i\\}_i`` denote the set of vertices of `X`.
Then any point ``p \\in \\mathbb{R}^n` of the convex polytope ``X`` is a convex
combination of its vertices, i.e., ``p = \\sum_{i} v_i α_i`` for some
(non-negative) coefficients ``\\{α_i\\}_i`` that add up to 1.
The algorithm chooses a random convex combination (the ``α_i``).
To produce such combination we apply the finite difference operator on a sorted
uniform sample over ``[0, 1]``; the method can be found in [1] and [2].

### Notes

The sampling is not uniform - points in the center of the polytope are more
likely to be sampled.

### References

[1] *Rubin, Donald B. The bayesian bootstrap. The annals of statistics (1981):
130-134.*

[2] https://cs.stackexchange.com/questions/3227/uniform-sampling-from-a-simplex/3229
"""
struct PolytopeSampler{S<:LazySet, D} <: Sampler
    X::S
end

PolytopeSampler(X::LazySet) = PolytopeSampler{typeof(X), DefaultUniform{eltype(X)}}(X)

set(sampler::PolytopeSampler) = sampler.X

# default sampler algorithms
_default_sampler(X::LazySet) = RejectionSampler
_default_sampler(X::LineSegment) = UniformSampler

"""
    _sample!(D::Vector{VN},
             sampler::RejectionSampler;
             rng::AbstractRNG=GLOBAL_RNG,
             seed::Union{Int, Nothing}=nothing) where {N, VN<:AbstractVector{N}}

Sample points using rejection sampling.

### Input

- `D`           -- output, vector of points
- `sampler`     -- Sampler from which the points are sampled
- `rng`         -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`        -- (optional, default: `nothing`) seed for reseeding

### Output

A vector of `num_samples` vectors.
"""
function _sample!(D::Vector{VN},
                  sampler::RejectionSampler;
                  rng::AbstractRNG=GLOBAL_RNG,
                  seed::Union{Int, Nothing}=nothing) where {N, VN<:AbstractVector{N}}
    rng = reseed(rng, seed)
    @inbounds for i in 1:length(D)
        w = rand.(Ref(rng), sampler.box_approx)
        while w ∉ sampler.X
            w = rand.(Ref(rng), sampler.box_approx)
        end
        D[i] = w
    end
    return D
end

"""
    _sample!(D::Vector{VN},
             sampler::UniformSampler{<:LineSegment, <:DefaultUniform};
             rng::AbstractRNG=GLOBAL_RNG,
             seed::Union{Int, Nothing}=nothing) where {N, VN<:AbstractVector{N}}

Sample points of a line segment using uniform sampling in-place.

### Input

- `D`           -- output, vector of points
- `sampler`     -- Sampler from which the points are sampled
- `rng`         -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`        -- (optional, default: `nothing`) seed for reseeding

### Output

A vector of `num_samples` vectors, where `num_samples` is the length of `D`.
"""
function _sample!(D::Vector{VN},
                  sampler::UniformSampler{<:LineSegment, <:DefaultUniform};
                  rng::AbstractRNG=GLOBAL_RNG,
                  seed::Union{Int, Nothing}=nothing) where {N, VN<:AbstractVector{N}}
    rng = reseed(rng, seed)
    U = DefaultUniform(zero(N), one(N))
    L = set(sampler)
    p = L.p
    q = L.q

    @inbounds for i in 1:length(D)
        λ = rand(rng, U)
        D[i] = p + λ * (q - p)
    end
    return D
end

function Base.rand(rng::AbstractRNG, U::DefaultUniform, n::Int)
    [rand(rng, U) for i in 1:n]
end

function rand!(x, rng::AbstractRNG, U::DefaultUniform)
    @inbounds for i in eachindex(x)
        x[i] = rand(rng, U)
    end
end

function _sample!(D::Vector{VN},
                  sampler::PolytopeSampler{<:LazySet, <:DefaultUniform};
                  rng::AbstractRNG=GLOBAL_RNG,
                  seed::Union{Int, Nothing}=nothing
                 ) where {N, VN<:AbstractVector{N}}
    rng = reseed(rng, seed)
    U = DefaultUniform(zero(N), one(N))
    P = set(sampler)
    vlist = vertices_list(P)
    m = length(vlist)

    # vector used to store the combination coefficients
    r = Vector{N}(undef, m-1)

    @inbounds for i in 1:length(D)
        # get a list of m uniform numbers (https://cs.stackexchange.com/a/3229)
        # and compute the corresponding linear combination in-place
        rand!(r, rng, U)
        sort!(r)
        D[i] = r[1] * vlist[1]  # r[1] - 0 == r[1]
        for j in 2:m-1
            α = r[j] - r[j-1]
            D[i] .+= α * vlist[j]
        end
        D[i] .+= (1 - r[m-1]) * vlist[m]
    end

    return D
end

# =============================
# Code requiring Distributions
# =============================

function load_distributions_samples()
return quote

using .Distributions: Uniform, Normal
import .Distributions

# ======================================================
# Sampling from a uniform distribution on balls/spheres
# ======================================================

"""
    _sample_unit_nsphere_muller!(D::Vector{Vector{N}}, n::Int, p::Int;
                                 rng::AbstractRNG=GLOBAL_RNG,
                                 seed::Union{Int, Nothing}=nothing) where {N}

Draw samples from a uniform distribution on an ``n``-dimensional unit
sphere using Muller's method.

### Input

- `D`    -- output, vector of points
- `n`    -- dimension of the sphere
- `p`    -- number of random samples
- `rng`  -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed` -- (optional, default: `nothing`) seed for reseeding

### Output

A vector of `nsamples` vectors.

### Algorithm

This function implements Muller's method of normalised Gaussians [1] to uniformly
sample over the ``n``-dimensional sphere ``S^n`` (which is the bounding surface
of the ``n``-dimensional unit ball).

Given ``n`` canonical Gaussian random variables ``Z₁, Z₂, …, Z_n``, the
distribution of the vectors

```math
\\dfrac{1}{α}\\left(z₁, z₂, …, z_n\\right)^T,
```
where ``α := \\sqrt{z₁² + z₂² + … + z_n²}``, is uniform over ``S^n``.

[1] Muller, Mervin E. *A note on a method for generating points uniformly on
    n-dimensional spheres.* Communications of the ACM 2.4 (1959): 19-20.
"""
function _sample_unit_nsphere_muller!(D::Vector{Vector{N}}, n::Int, p::Int;
                                      rng::AbstractRNG=GLOBAL_RNG,
                                      seed::Union{Int, Nothing}=nothing) where {N}
    rng = reseed(rng, seed)
    Zdims = [Normal() for _ in 1:n] # normal distributions, one for each dimension
    v = Vector{N}(undef, n) # sample direction
    @inbounds for j in 1:p
        α = zero(N)
        for i in 1:n
            v[i] = rand(rng, Zdims[i])
            α += v[i]^2
        end
        D[j] = v ./ sqrt(α)
    end
    return D
end

"""
    _sample_unit_nball_muller!(D::Vector{Vector{N}}, n::Int, p::Int;
                               rng::AbstractRNG=GLOBAL_RNG,
                               seed::Union{Int, Nothing}=nothing) where {N}

Draw samples from a uniform distribution on an ``n``-dimensional unit ball
using Muller's method.

### Input

- `D`    -- output, vector of points
- `n`    -- dimension of the ball
- `p`    -- number of random samples
- `rng`  -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed` -- (optional, default: `nothing`) seed for reseeding

### Output

A vector of `nsamples` vectors.

### Algorithm

This function implements Muller's method of normalised Gaussians [1] to
sample from the interior of the ball.

Given ``n`` Gaussian random variables ``Z₁, Z₂, …, Z_n`` and a uniformly
distributed random variable ``r`` with support in ``[0, 1]``, the distribution
of the vectors

```math
\\dfrac{r^{1/n}}{α} \\left(z₁, z₂, …, z_n\\right)^T,
```
where ``α := \\sqrt{z₁² + z₂² + … + z_n²}``, is uniform over the
``n``-dimensional unit ball.

[1] Muller, Mervin E. *A note on a method for generating points uniformly on
    n-dimensional spheres.* Communications of the ACM 2.4 (1959): 19-20.
"""
function _sample_unit_nball_muller!(D::Vector{Vector{N}}, n::Int, p::Int;
                                    rng::AbstractRNG=GLOBAL_RNG,
                                    seed::Union{Int, Nothing}=nothing) where {N}

    rng = reseed(rng, seed)
    Zdims = [Normal() for _ in 1:n] # normal distributions, one for each dimension
    Zrad = Uniform() # distribution to pick random radius
    one_over_n = one(N)/n
    v = Vector{N}(undef, n) # sample direction
    @inbounds for j in 1:p
        α = zero(N)
        for i in 1:n
            v[i] = rand(rng, Zdims[i])
            α += v[i]^2
        end
        r = rand(rng, Zrad)
        β = r^one_over_n / sqrt(α)
        D[j] = v .* β
    end
    return D
end

end # quote
end # function load_distributions_samples()
