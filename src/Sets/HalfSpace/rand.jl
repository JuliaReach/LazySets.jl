"""
    rand(::Type{HalfSpace}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

Create a random half-space.

### Input

- `HalfSpace` -- type for dispatch
- `N`         -- (optional, default: `Float64`) numeric type
- `dim`       -- (optional, default: 2) dimension
- `rng`       -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`      -- (optional, default: `nothing`) seed for reseeding

### Output

A random half-space.

### Algorithm

All numbers are normally distributed with mean 0 and standard deviation 1.
Additionally, the constraint `a` is nonzero.
"""
function rand(::Type{HalfSpace};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing)
    rng = reseed!(rng, seed)
    a = randn(rng, N, dim)
    while iszero(a)
        a = randn(rng, N, dim)  # COV_EXCL_LINE
    end  # COV_EXCL_LINE
    b = randn(rng, N)
    return HalfSpace(a, b)
end
