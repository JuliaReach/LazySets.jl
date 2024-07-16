"""
    rand(::Type{Star}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

Create a random star.

### Input

- `Star` -- type for dispatch
- `N`    -- (optional, default: `Float64`) numeric type
- `dim`  -- (optional, default: 2) dimension
- `rng`  -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed` -- (optional, default: `nothing`) seed for reseeding
- `P`    -- (optional, default: a random `HPolytope`) predicate

### Output

A random star.

### Algorithm

By default we generate a random `HPolytope` of dimension `dim` as predicate.
Alternatively the predicate can be passed.

All numbers are normally distributed with mean 0 and standard deviation 1.
"""
function rand(::Type{Star};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing,
              P::AbstractPolyhedron=rand(HPolytope; N=N, dim=dim, rng=rng, seed=seed))
    rng = reseed!(rng, seed)
    c = randn(rng, N, dim)
    V = randn(rng, N, dim, LazySets.dim(P))
    return Star(c, V, P)
end
