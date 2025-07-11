"""
# Extended help

    rand(::Type{MatrixZonotope}; [N]::Type{<:Real}=Float64, [dim]::Tuple{Int,Int}=(2, 2),
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing,
         [num_generators]::Int=-1)

### Algorithm

All numbers are normally distributed with mean 0 and standard deviation 1.

The number of generators can be controlled with the argument `num_generators`.
For a negative value we choose a random number in the range `1:maximum(dim)`.
"""
function rand(::Type{MatrixZonotope};
              N::Type{<:Real}=Float64,
              dim::Tuple{Int,Int}=(2, 2),
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing,
              num_generators::Int=-1)
    rng = reseed!(rng, seed)
    center = randn(rng, N, dim)
    if num_generators < 0
        num_generators = rand(rng, 1:maximum(dim))
    end

    generators = Vector{Matrix{N}}(undef, num_generators)
    @inbounds for i in 1:num_generators
        generators[i] = randn(rng, N, dim)
    end
    return MatrixZonotope(center, generators)
end
