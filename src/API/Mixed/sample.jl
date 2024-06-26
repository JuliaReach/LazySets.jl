"""
    sample(X::LazySet, [m]::Int=1;
           [rng]::AbstractRNG=GLOBAL_RNG,
           [seed]::Union{Int,Nothing}=nothing)

Compute random samples from a set.

### Input

- `X`    -- set
- `m`    -- (optional; default: 1) number of random samples
- `rng`  -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed` -- (optional, default: `nothing`) seed for reseeding

### Output

A vector of `m` elements in `X` if `X` is nonempty, and an error otherwise.
"""
function sample(::LazySet, ::Int=1) end
