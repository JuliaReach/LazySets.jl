module LazySetsDistributionsExt

using Distributions: Normal, UnivariateDistribution
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Distribution: reseed!
import LazySets: RejectionSampler, _sample_unit_nsphere_muller!

include("Distributions/Ball2Ext.jl")

function RejectionSampler(distr::UnivariateDistribution; tight::Bool=false)
    return RejectionSampler([distr]; tight=tight)
end

"""
    _sample_unit_nsphere_muller!(D::Vector{Vector{N}}, n::Int, p::Int;
                                 [rng]::AbstractRNG=GLOBAL_RNG,
                                 [seed]::Union{Int, Nothing}=nothing) where {N}

Draw samples from a uniform distribution on an ``n``-dimensional unit sphere
using Muller's method.

### Input

- `D`    -- output, vector of points
- `n`    -- dimension of the sphere
- `p`    -- number of random samples
- `rng`  -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed` -- (optional, default: `nothing`) seed for reseeding

### Output

The modified vector `D`.

### Algorithm

This function implements Muller's method of normalized Gaussians [Muller59](@citet)
to uniformly sample over the ``n``-dimensional sphere ``S^n`` (which is the
bounding surface of the ``n``-dimensional unit ball).

Given ``n`` canonical Gaussian random variables ``Z₁, Z₂, …, Z_n``, the
distribution of the vectors

```math
\\dfrac{1}{α}\\left(z₁, z₂, …, z_n\\right)^T,
```
where ``α := \\sqrt{z₁² + z₂² + … + z_n²}``, is uniform over ``S^n``.
"""
function _sample_unit_nsphere_muller!(D::Vector{Vector{N}}, n::Int, p::Int;
                                      rng::AbstractRNG=GLOBAL_RNG,
                                      seed::Union{Int,Nothing}=nothing) where {N}
    rng = reseed!(rng, seed)
    Zdims = [Normal() for _ in 1:n] # normal distributions for each dimension
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

end  # module
