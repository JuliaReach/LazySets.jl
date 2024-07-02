"""
    rand(::Type{Ball1}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing
        )

Create a random ball in the 1-norm.

### Input

- `Ball1` -- type for dispatch
- `N`     -- (optional, default: `Float64`) numeric type
- `dim`   -- (optional, default: 2) dimension
- `rng`   -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`  -- (optional, default: `nothing`) seed for reseeding

### Output

A random ball in the 1-norm.

### Algorithm

All numbers are normally distributed with mean 0 and standard deviation 1.
Additionally, the radius is nonnegative.
"""
function rand(::Type{Ball1};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing)
    rng = reseed!(rng, seed)
    center = randn(rng, N, dim)
    radius = abs(randn(rng, N))
    return Ball1(center, radius)
end
