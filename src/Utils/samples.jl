# ========================
# Sampling from a LazySet
# ========================

"""
    AbstractSampler

Abstract type for defining new sampling methods.

### Notes

All subtypes should implement a `sample!(D, X, ::Method)` method where the
first argument is the output (vector of vectors), the second argument is the
set to be sampled, and the third argument is the sampler instance.
"""
abstract type AbstractSampler end

"""
    sample(X::LazySet{N}, num_samples::Int;
           [sampler]=_default_sampler(X),
           [rng]::AbstractRNG=GLOBAL_RNG,
           [seed]::Union{Int, Nothing}=nothing,
           [include_vertices]=false,
           [VN]=Vector{N}) where {N}

Random sampling of an arbitrary set `X`.

### Input

- `X`           -- set to be sampled
- `num_samples` -- number of random samples
- `sampler`     -- (optional, default: `_default_sampler(X)`) the sampler used;
                   falls back to `RejectionSampler`
- `rng`         -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`        -- (optional, default: `nothing`) seed for reseeding
- `include_vertices` -- (optional, default: `false`) option to include the vertices of `X`
- `VN`          -- (optional, default: `Vector{N}`) vector type of the sampled points

### Output

A vector of `num_samples` vectors.
If `num_samples` is not passed, the result is just one sample (not wrapped in a
vector).

### Algorithm

See the documentation of the respective `Sampler`.

### Notes

If `include_vertices == true`, we include all vertices computed with `vertices`.
Alternatively if a number ``k`` is passed, we plot the first ``k`` vertices
returned by `vertices(X)`.
"""
function sample(X::LazySet{N}, num_samples::Int;
                sampler=_default_sampler(X),
                rng::AbstractRNG=GLOBAL_RNG,
                seed::Union{Int, Nothing}=nothing,
                include_vertices=false,
                VN=Vector{N}) where {N}

    D = Vector{VN}(undef, num_samples) # preallocate output
    sample!(D, X, sampler; rng=rng, seed=seed)

    if include_vertices != false
        k = (include_vertices isa Bool) ? Inf : include_vertices
        for v in vertices(X)
            push!(D, v)
            k -= 1
            if k <= 0
                break
            end
        end
    end

    return D
end

# without argument, returns a single element (instead of a singleton)
function sample(X::LazySet; kwargs...)
    return sample(X, 1; kwargs...)[1]
end

# =====================
# Uniform distribution
# =====================

# represents a uniform distribution over the interval [a, b]
# using `rand` from the Julia standard library
struct DefaultUniform{N}
    a::N
    b::N
end

function Base.rand(rng::AbstractRNG, U::DefaultUniform)
    r = rand(rng)
    Δ = U.b - U.a
    return Δ * r + U.a
end

function Base.rand(rng::AbstractRNG, U::DefaultUniform, n::Int)
    return [rand(rng, U) for i in 1:n]
end

function Base.rand(rng::AbstractRNG, U::AbstractVector{<:DefaultUniform})
    return rand.(Ref(rng), U)
end

function rand!(x, rng::AbstractRNG, U::DefaultUniform)
    @inbounds for i in eachindex(x)
        x[i] = rand(rng, U)
    end
    return x
end

# ===================
# Rejection Sampling
# ===================

"""
    RejectionSampler{D} <: Sampler

Type used for rejection sampling of an arbitrary set `X`.

### Fields

- `distribution` -- (optional, default: `DefaultUniform`) distribution from which
                    the sample is drawn
- `tight`        -- (optional, default: `false`) set to `true` if the support of
                    the distribution is known to coincide with the set `X`

### Algorithm

Draw a sample ``x`` from a given distribution of a box-overapproximation of the
original set ``X`` in all ``n`` dimensions. The function rejects a drawn sample
``x`` and redraws as long as the sample is not contained in the original set
``X``, i.e., ``x ∉ X``.
"""
struct RejectionSampler{D} <: AbstractSampler
    distribution::D
    tight::Bool
end

function RejectionSampler(distr; tight::Bool=false)
    return RejectionSampler(distr, tight)
end

function RejectionSampler(distr::DefaultUniform; tight::Bool=false)
    return RejectionSampler([distr], tight)
end

function RejectionSampler(X::LazySet, distribution=DefaultUniform; tight::Bool=false)
    # define the support of the distribution as the smallest box enclosing X
    n = dim(X)
    B = box_approximation(X)

    # distribution over B
    distr = [distribution(low(B, i), high(B, i)) for i in 1:n]

    return RejectionSampler(distr, tight)
end

# the support of this distribution is always tight wrt X
function RejectionSampler(X::AbstractHyperrectangle)
    n = dim(X)
    distr = [DefaultUniform(low(X, i), high(X, i)) for i in 1:n]
    return RejectionSampler(distr, true)
end

# ambiguity fix
function RejectionSampler(X::LazySet, tight::Bool)
    RejectionSampler(X, DefaultUniform; tight=tight)
end
# default sampling for LazySets
_default_sampler(X::LazySet) = RejectionSampler(X)
_default_sampler(X::LineSegment{N}) where {N} = RejectionSampler(DefaultUniform(zero(N), one(N)), true)

function sample!(D::Vector{VN}, X::LazySet, sampler::RejectionSampler;
                 rng::AbstractRNG=GLOBAL_RNG,
                 seed::Union{Int, Nothing}=nothing) where {N, VN<:AbstractVector{N}}

    U = sampler.distribution
    rng = reseed(rng, seed)
    @inbounds for i in 1:length(D)
        w = rand(rng, U)

        if !(sampler.tight)
            while w ∉ X
                w = rand(rng, U)
            end
        end
        D[i] = w
    end
    return D
end

function sample!(D::Vector{VN}, L::LineSegment, sampler::RejectionSampler{<:DefaultUniform};
                 rng::AbstractRNG=GLOBAL_RNG,
                 seed::Union{Int, Nothing}=nothing) where {N, VN<:AbstractVector{N}}

    rng = reseed(rng, seed)
    U = sampler.distribution
    @assert U.a >= zero(N) && U.b <= one(N)
    p = L.p
    q = L.q

    @inbounds for i in 1:length(D)
        λ = rand(rng, U)
        D[i] = p + λ * (q - p)
    end
    return D
end

"""
    RandomWalkSampler{D} <: Sampler

Type used for sampling of a convex polytope `X` using its vertex representation.

### Fields

- `X` -- convex polytope to be sampled

### Algorithm

Choose a random convex combination of the vertices of `X`.

Let ``V = \\{v_i\\}_i`` denote the set of vertices of `X`.
Then any point ``p \\in \\mathbb{R}^n`` of the convex polytope ``X`` is a convex
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
struct RandomWalkSampler <: AbstractSampler
    #
end

function sample!(D::Vector{VN}, X::LazySet, sampler::RandomWalkSampler;
                 rng::AbstractRNG=GLOBAL_RNG,
                 seed::Union{Int, Nothing}=nothing) where {N, VN<:AbstractVector{N}}
    rng = reseed(rng, seed)
    U = DefaultUniform(zero(N), one(N))
    vlist = vertices_list(X)
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

using .Distributions: Uniform, Normal, Distribution, UnivariateDistribution
import .Distributions

RejectionSampler(distr::UnivariateDistribution; tight::Bool=false) = RejectionSampler([distr], tight=tight)

function Base.rand(rng::AbstractRNG, U::AbstractVector{<:UnivariateDistribution})
    return rand.(Ref(rng), U)
end

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
