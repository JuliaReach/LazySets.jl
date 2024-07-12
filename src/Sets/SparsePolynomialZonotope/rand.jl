"""
    rand(::Type{SparsePolynomialZonotope}; [N]::Type{<:Real}=Float64,
         [dim]::Int=2, [nparams]::Int=2, [maxdeg]::Int=3,
         [num_dependent_generators]::Int=-1,
         [num_independent_generators]::Int=-1, [rng]::AbstractRNG=GLOBAL_RNG,
         [seed]::Union{Int, Nothing}=nothing)

Create a random sparse polynomial zonotope.

### Input

- `SparsePolynomialZonotope`   -- type for dispatch
- `N`                          -- (optional, default: `Float64`) numeric type
- `dim`                        -- (optional, default: 2) dimension
- `nparams`                    -- (optional, default: 2) number of parameters
- `maxdeg`                     -- (optional, default: 3) maximum degree for each
                                  parameter
- `num_dependent_generators`   -- (optional, default: `-1`) number of dependent
                                  generators (see comment below)
- `num_independent_generators` -- (optional, default: `-1`) number of
                                  independent generators (see comment below)
- `rng`                        -- (optional, default: `GLOBAL_RNG`) random
                                  number generator
- `seed`                       -- (optional, default: `nothing`) seed for
                                  reseeding

### Output

A random sparse polynomial zonotope.

### Algorithm

All numbers are normally distributed with mean 0 and standard deviation 1.

The number of generators can be controlled with the arguments
`num_dependent_generators` and `num_dependent_generators`.
For a negative value we choose a random number in the range `dim:2*dim` (except
if `dim == 1`, in which case we only create a single generator). Note that the
final number of generators may be lower if redundant monomials are generated.
"""
function rand(::Type{SparsePolynomialZonotope};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              nparams::Int=2,
              maxdeg::Int=3,
              num_dependent_generators::Int=-1,
              num_independent_generators::Int=-1,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing)
    rng = reseed!(rng, seed)

    if num_independent_generators < 0
        num_independent_generators = (dim == 1) ? 1 : rand(rng, dim:(2 * dim))
    end
    GI = randn(rng, N, dim, num_independent_generators)

    SSPZ = rand(SimpleSparsePolynomialZonotope; N=N, dim=dim, nparams=nparams,
                maxdeg=maxdeg, rng=rng, num_generators=num_dependent_generators)

    return SparsePolynomialZonotope(center(SSPZ), genmat(SSPZ), GI, expmat(SSPZ))
end
