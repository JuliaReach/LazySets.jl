"""
    rand(::Type{Singleton}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

Create a random singleton.

### Input

- `Singleton` -- type for dispatch
- `N`         -- (optional, default: `Float64`) numeric type
- `dim`       -- (optional, default: 2) dimension
- `rng`       -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`      -- (optional, default: `nothing`) seed for reseeding

### Output

A random singleton.

### Algorithm

The element is a normally distributed vector with entries of mean 0 and standard
deviation 1.
"""
function rand(::Type{Singleton};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing)
    rng = reseed!(rng, seed)
    element = randn(rng, N, dim)
    return Singleton(element)
end
