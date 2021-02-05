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
    sample(X::LazySet{N}, num_samples::Int, [sampler]::Type{Sampler}=_default_sampler(X);
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
function sample(X::LazySet{N}, num_samples::Int, sampler::Type{<:Sampler}=_default_sampler(X);
                rng::AbstractRNG=GLOBAL_RNG,
                seed::Union{Int, Nothing}=nothing,
                VN=Vector{N}) where {N}

    S = sampler(X)
    return sample(X, num_samples, S, rng=rng, seed=seed, VN=VN)
end

function sample(X::LazySet{N}, num_samples::Int, sampler::Sampler;
                rng::AbstractRNG=GLOBAL_RNG,
                seed::Union{Int, Nothing}=nothing,
                VN=Vector{N}) where {N}

    D = Vector{VN}(undef, num_samples) # preallocate output
    _sample!(D, sampler; rng=rng, seed=seed)
    return D
end

# without argument, returns a single sample (instead of a vector containing a single vector)
sample(X::LazySet; kwargs...) = sample(X, 1; kwargs...)[1]
sample(X::LazySet, sampler::Union{Type{Sampler}, Sampler}; kwargs...) = sample(X, 1, sampler; kwargs...)[1]

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


# =========================
# Hit and Run Monte Carlo
# =========================

struct HitAndRun{S, D, PT, L} <: Sampler
    X::S      #  set
    p::PT     # starting point
    thin::L   # thin level
end

set(sampler::HitAndRun) = sampler.X

# FIXME an_element(X) gives a point in the boundary, but the algorithm works
# better with an interior point
function HitAndRun(X::LazySet; p=an_element(X), thin=10)
    return HitAndRun{typeof(X), DefaultUniform{eltype(X)}, typeof(p), typeof(thin)}(X, p, thin)
end

# compute a unitay random direction in the given dimensions
function _unitary_random_direction(n)
    return normalize!(randn(n))
end

# distance we have to travel in direction dir,
# from point p, to reach a given hyperplane a^T x = b
function _distance(a, b, p, dir)
    α = dot(a, dir)
    LazySets.isapproxzero(α) && throw(ArgumentError("the direction is parallel to the hyperplane"))
    λ = (b - dot(a, p)) / α
    return λ
end

function _distance(H::Hyperplane, p, dir)
    return _distance(H.a, H.b, p, dir)
end

function _sample!(D::Vector{VN},
                  sampler::HitAndRun;
                  rng::AbstractRNG=GLOBAL_RNG,
                  seed::Union{Int, Nothing}=nothing) where {N, VN<:AbstractVector{N}}

    # initialization
    rng = reseed(rng, seed)
    X = set(sampler)
    p = an_element(X)
    n = dim(X)
    clist = constraints_list(X)
    D[1] = p

    thin = sampler.thin
    if thin == 1
        for i in 2:length(D)
            D[i] = hitandrun(clist, D[i-1])
        end

    else
        aux = similar(p)
        for i in 2:length(D)
            D[i] = similar(p)
            copy!(aux, D[i-1])
            for _ in 1:thin
                hitandrun!(D[i], clist, aux)
                copy!(aux, D[i])
            end
        end
    end
    return D
end

hitandrun(clist::Vector{<:HalfSpace}, p) = hitandrun!(similar(p), clist, p)

function hitandrun!(q, clist::Vector{HT}, p; MAXITER=1000) where {N, VT, HT<:HalfSpace{N, VT}}

    n = dim(first(clist))
    m = length(clist)
    λ = Vector{N}(undef, m)
    k = 1
    while k <= MAXITER
        dir = _unitary_random_direction(n)
        @inbounds for (i, c) in enumerate(clist)
            λ[i] = _distance(c.a, c.b, p, dir)
        end

        λ₋ = -Inf   # maximum of negative coeffs
        λ₊ = Inf # minimum of positive coeffs
        for λi in λ
            if λi < zero(N) && λi > λ₋
                λ₋ = λi
            end
            if λi > zero(N) && λi < λ₊
                λ₊ = λi
            end
        end

        # build pointf
        U = DefaultUniform(λ₋, λ₊)
        μ = rand(U)
        q .= p + μ * dir

        # check membership
        found = true
        for c in clist
            if !_leq(dot(c.a, q), c.b)
                found = false
                break
            end
        end
        found && break
        k += 1
    end
    if k > MAXITER
        @warn "maximum number of hit and run iterations reached, try increasing `MAXITER` or choose another starting point"
    end
    return q
end
