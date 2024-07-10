"""
    rand(::Type{Universe}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

Create a universe (note that there is nothing to randomize).

### Input

- `Universe` -- type for dispatch
- `N`        -- (optional, default: `Float64`) numeric type
- `dim`      -- (optional, default: 2) dimension
- `rng`      -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`     -- (optional, default: `nothing`) seed for reseeding

### Output

The (only) universe of the given numeric type and dimension.
"""
function rand(::Type{Universe};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing)
    rng = reseed!(rng, seed)
    return Universe{N}(dim)
end
