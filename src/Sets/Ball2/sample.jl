"""
# Extended help

    sample(B::Ball2, [nsamples]::Int;
           [rng]::AbstractRNG=GLOBAL_RNG,
           [seed]::Union{Int, Nothing}=nothing)

### Algorithm

Random sampling with uniform distribution in `B` is computed using Muller's method
of normalized Gaussians. This method requires the package `Distributions`.
See [`_sample_unit_nball_muller!`](@ref) for implementation details.
"""
function sample(B::Ball2, nsamples::Int;
                rng::AbstractRNG=GLOBAL_RNG,
                seed::Union{Int,Nothing}=nothing)
    n = dim(B)
    N = eltype(B)
    D = Vector{Vector{N}}(undef, nsamples) # preallocate output
    _sample_unit_nball_muller!(D, n, nsamples; rng=rng, seed=seed)

    # customize for the given ball
    r, c = B.radius, B.center
    @inbounds for i in 1:nsamples
        axpby!(one(N), c, r, D[i])
    end
    return D
end

# see ext/Distributions/DistributionsBall2Ext.jl
function _sample_unit_nball_muller!(D, n, p; rng=GLOBAL_RNG, seed=nothing)
    mod = Base.get_extension(@__MODULE__, :DistributionsExt)
    require(mod, :Distributions; fun_name="_sample_unit_nball_muller!")
    error()
end
