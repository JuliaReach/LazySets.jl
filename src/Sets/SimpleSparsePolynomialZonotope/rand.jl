"""
    rand(::Type{SimpleSparsePolynomialZonotope};
         [N]::Type{<:Real}=Float64, [dim]::Int=2, [nparams]::Int=2,
         [maxdeg]::Int=3, [num_generators]::Int=-1,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

Create a random simple sparse polynomial zonotope.

### Input

- `Zonotope`       -- type for dispatch
- `N`              -- (optional, default: `Float64`) numeric type
- `dim`            -- (optional, default: 2) dimension
- `nparams`        -- (optional, default: 2) number of parameters
- `maxdeg`         -- (optional, default: 3) maximum degree for each parameter
- `rng`            -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`           -- (optional, default: `nothing`) seed for reseeding
- `num_generators` -- (optional, default: `-1`) number of generators of the
                      zonotope (see comment below)

### Output

A random simple sparse polynomial zonotope.

### Algorithm

All numbers are normally distributed with mean 0 and standard deviation 1.

The number of generators can be controlled with the argument `num_generators`.
For a negative value we choose a random number in the range `dim:2*dim` (except
if `dim == 1`, in which case we only create a single generator). Note that the
final number of generators may be lower if redundant monomials are generated.
"""
function rand(::Type{SimpleSparsePolynomialZonotope};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              nparams::Int=2,
              maxdeg::Int=3,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing,
              num_generators::Int=-1)
    rng = reseed!(rng, seed)
    center = randn(rng, N, dim)
    if num_generators < 0
        num_generators = (dim == 1) ? 1 : rand(rng, dim:(2 * dim))
    end
    generators = randn(rng, N, dim, num_generators)
    expmat = rand(rng, 0:maxdeg, nparams, num_generators)
    SSPZ = SimpleSparsePolynomialZonotope(center, generators, expmat)
    return remove_redundant_generators(SSPZ)
end
