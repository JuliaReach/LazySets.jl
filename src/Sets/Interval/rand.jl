"""
    rand(::Type{Interval}; [N]::Type{<:Real}=Float64, [dim]::Int=1,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

Create a random interval.

### Input

- `Interval` -- type for dispatch
- `N`        -- (optional, default: `Float64`) numeric type
- `dim`      -- (optional, default: 1) dimension
- `rng`      -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`     -- (optional, default: `nothing`) seed for reseeding

### Output

A random interval.

### Algorithm

All numbers are normally distributed with mean 0 and standard deviation 1.
"""
function rand(::Type{Interval};
              N::Type{<:Real}=Float64,
              dim::Int=1,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing)
    @assert dim == 1 "cannot create a random Interval of dimension $dim"
    rng = reseed!(rng, seed)
    x = randn(rng, N)
    y = randn(rng, N)
    return x < y ? Interval(x, y) : Interval(y, x)
end
