"""
    sample(X::LazySet, m::Int;
           [rng]::AbstractRNG=GLOBAL_RNG,
           [seed]::Union{Int,Nothing}=nothing)

Compute samples from a set.

### Input

- `X`    -- set
- `m`    -- number of random samples
- `rng`  -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed` -- (optional, default: `nothing`) seed for reseeding

### Output

A vector of `m` elements in `X`.
"""
function sample(::LazySet, ::Int) end
