function load_distributions_samples()
return quote

using .Distributions: Distribution, Uniform, Normal
import .Distributions

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

## Rejection Sampling

"""
    sample(X::LazySet{N}, num_samples::Int;
           rng::AbstractRNG=GLOBAL_RNG,
           seed::Union{Int, Nothing}=nothing) where {N}

Sampling of an arbitrary `LazySet` `X`.

### Input

- `X`           -- lazyset
- `num_samples` -- number of random samples
- `sampler`     -- Sampler used (default: `RejectionSampler`)
- `rng`         -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`        -- (optional, default: `nothing`) seed for reseeding

### Output

A vector of `num_samples` vectors.
If `num_samples` is not passed one sample as a single vector is returned.

### Algorithm

`RejectionSampler`:
Draw a sample `x` from a uniform distribution of a box-overapproximation of the
original set `X` in all `n` dimensions. The function rejects a drawn sample `x`
and redraws as long as the sample is not contained in the original set `X`,
i.e., `x ∉ X`.
"""
function sample(X::LazySet{N}, num_samples::Int;
                sampler=RejectionSampler,
                rng::AbstractRNG=GLOBAL_RNG,
                seed::Union{Int, Nothing}=nothing) where {N<:Real}
    @assert isbounded(X) "this function requires that the set `X` is bounded"*
                            ", but it is not"

    D = Vector{Vector{N}}(undef, num_samples) # preallocate output
    _sample_from_sampler!(D, sampler(X); rng=rng, seed=seed)
    return D
end

function sample(X::LazySet{N}; kwargs...) where {N<:Real}
    return sample(X, 1; kwargs...)[1]
end

"""
    Sampler

Abstract type for defining new sample methods.

### Notes

All subtypes should implement a `_sample_from_sampler!` method.
"""
abstract type Sampler end

"""
    RejectionSampler{S<:LazySet, D<:Distribution} <: Sampler

Type used for rejection sampling of an arbitrary `LazySet` `X`.

### Fields

- `X`           -- lazyset
- `box_approx`  -- Distribution from which the sample is drawn
"""
struct RejectionSampler{S<:LazySet, D<:Distribution} <:Sampler
    X::S
    box_approx::Vector{D}
end

function RejectionSampler(X, distribution=Uniform)
    B = box_approximation(X)
    canonical_support =  hcat(low(B), high(B))
    dims = size(canonical_support,1)
    box_approx = [distribution(canonical_support[i,:]...) for i = 1:dims]
    return RejectionSampler(X, box_approx)
end


"""
    _sample_from_sampler!(D::Vector{Vector{N}},
           sampler::RejectionSampler;
           rng::AbstractRNG=GLOBAL_RNG,
           seed::Union{Int, Nothing}=nothing) where {N<:Real}

Sampling points from a `Sampler`.

### Input

- `D`           -- output, vector of points
- `sampler`     -- Sampler from which the points are sampled
- `rng`         -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`        -- (optional, default: `nothing`) seed for reseeding

### Output

A vector of `num_samples` vectors.
"""
function _sample_from_sampler!(D::Vector{Vector{N}},
               sampler::RejectionSampler;
               rng::AbstractRNG=GLOBAL_RNG,
               seed::Union{Int, Nothing}=nothing) where {N<:Real}
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

end # quote
end # function load_distributions_samples()
