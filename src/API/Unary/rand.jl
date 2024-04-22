"""
    rand(T::Type{<:LazySet}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing
        )

Create a random set of the given set type.

### Input

- `T`    -- set type
- `N`    -- (optional, default: `Float64`) numeric type
- `dim`  -- (optional, default: 2) dimension
- `rng`  -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed` -- (optional, default: `nothing`) seed for reseeding

### Output

A random set of the given set type.
"""
function rand(::Type{<:LazySet}) end
