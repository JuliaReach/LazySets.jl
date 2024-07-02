"""
    sample(B::Ball2{N}, [nsamples]::Int;
           [rng]::AbstractRNG=GLOBAL_RNG,
           [seed]::Union{Int, Nothing}=nothing) where {N}

Return samples from a uniform distribution on the given ball in the 2-norm.

### Input

- `B`        -- ball in the 2-norm
- `nsamples` -- number of random samples
- `rng`      -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`     -- (optional, default: `nothing`) seed for reseeding

### Output

A linear array of `nsamples` elements drawn from a uniform distribution in `B`.

### Algorithm

Random sampling with uniform distribution in `B` is computed using Muller's method
of normalized Gaussians. This function requires the package `Distributions`.
See `_sample_unit_nball_muller!` for implementation details.
"""
function sample(B::Ball2{N}, nsamples::Int;
                rng::AbstractRNG=GLOBAL_RNG,
                seed::Union{Int,Nothing}=nothing) where {N}
    require(@__MODULE__, :Distributions; fun_name="sample")
    n = dim(B)
    D = Vector{Vector{N}}(undef, nsamples) # preallocate output
    _sample_unit_nball_muller!(D, n, nsamples; rng=rng, seed=seed)

    # customize for the given ball
    r, c = B.radius, B.center
    @inbounds for i in 1:nsamples
        axpby!(one(N), c, r, D[i])
    end
    return D
end
