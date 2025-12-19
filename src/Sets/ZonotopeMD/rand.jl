"""
# Extended help

    rand(::Type{ZonotopeMD}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing,
         [num_generators]::Int=-1)

### Algorithm

All numbers are normally distributed with mean 0 and standard deviation 1.

The number of non-axis-aligned generators can be controlled with the argument
`num_generators`. For a negative value, we choose a random number in the range
`dim:2*dim` (except if `dim == 1`, in which case we only create a single
generator).
"""
function rand(::Type{ZonotopeMD};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing,
              num_generators::Int=-1)
    rng = reseed!(rng, seed)
    center = randn(rng, N, dim)
    if num_generators < 0
        num_generators = (dim == 1) ? 0 : rand(rng, dim:(2 * dim))
    end
    M = randn(rng, N, dim, num_generators)
    D = randn(rng, N, dim)
    return ZonotopeMD(center, M, D)
end
