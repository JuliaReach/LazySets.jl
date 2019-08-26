function load_distributions_samples()
return quote

using .Distributions: Distribution, Uniform, Normal

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
    _canonical_length(X::LazySet)

Returns the support function of the given set along the positive and negative canonical directions.

### Inputs

- `X` - convex set

### Outputs

A matrix with `n = dims(X)` rows and two columns. Each row stands for
one dimension of `X` whereas the first column is the minimum and the second
column is the maximum value of the corresponding dimension.
"""
function _canonical_length(X::LazySet{N}) where {N<:Real}
    dims = dim(X)
    x = Matrix{N}(undef, dims, 2)
    for j = 1:dims
        ej = SingleEntryVector(j, dims, one(N))
        x[j,:] = [-ρ(-ej, X), ρ(ej, X)]
    end
    return x
end

"""
    sample(X::LazySet{N}, num_samples::Int;
           rng::AbstractRNG=GLOBAL_RNG,
           seed::Union{Int, Nothing}=nothing) where {N}

Rejection sampling of an arbitrary LazySet `X` for which the support value
function is defined.

Draw a sample `x` from a uniform distribution of a box-overapproximation of the
original set `X` in all `n` dimensions. The function rejects a drawn sample `x`
and redraws as long as the sample is not contained in the original set `X`,
i.e., `x ∉ X`.

### Input

- `X`           -- lazyset
- `num_samples`    -- number of random samples
- `rng`         -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`        -- (optional, default: `nothing`) seed for reseeding

### Output

A vector of `num_samples` vectors.
"""
function Base.sample(X::LazySet{N}, num_samples::Int;
                rng::AbstractRNG=GLOBAL_RNG,
                seed::Union{Int, Nothing}=nothing) where {N<:Real}
    @assert isbounded(X) "this function requires that the set `X` is bounded, but it is not"

    D = Vector{Vector{N}}(undef, num_samples) # preallocate output
    sample_from_set!(D, Sampler(X); rng=rng, seed=seed)
    return D
end

function Base.sample(X::LazySet{N}; kwargs...) where {N<:Real}
    return sample(X, 1; kwargs...)[1]
end

struct Sampler{S<:LazySet, D<:Distribution}
    X::S
    box_approx::Vector{D}
end

function Sampler(X, distribution=Uniform)
    box = _canonical_length(X)
    sampler = [distribution(box[i,:]...) for i = 1:size(box,1)]
    return Sampler(X, sampler)
end

function sample_from_set!(D::Vector{Vector{N}},
               sampler::Sampler;
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
    nothing
end

end # quote
end # function load_distributions_samples()
