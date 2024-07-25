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
- `P`    -- (optional, default: `nothing`) predicate

### Output

A random star.

### Algorithm

A predicate `P` can be passed directly. If `P` is `nothing` (default), we
generate a random `HPolytope` of dimension `dim`.

All numbers are normally distributed with mean 0 and standard deviation 1.
"""
function rand(::Type{Star};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing,
              P::Union{AbstractPolyhedron,Nothing}=nothing)
    rng = reseed!(rng, seed)
    c = randn(rng, N, dim)
    V = randn(rng, N, dim, LazySets.dim(P))
    if isnothing(P)
        P = rand(HPolytope; N=N, dim=dim, rng=rng, seed=seed)
    end
    return Star(c, V, P)
end
