import LazySets

using LazySets: center, expmat, genmat, genmat_dep, genmat_indep, indexvector,
                ngens_dep, ngens_indep, reduce_order
using LazySets.Approximations: overapproximate
using LazySets.SimpleSparsePolynomialZonotopeModule: SimpleSparsePolynomialZonotope
using LazySets.SparsePolynomialZonotopeModule: SparsePolynomialZonotope
using LazySets.ZonotopeModule: Zonotope
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Distribution: reseed!
import Base: rand
import LazySets.SparsePolynomialZonotopeModule: _absorb_generators_spz

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

# absorb all generators whose norm is below a given threshold with a box
function _absorb_generators_spz(P::SparsePolynomialZonotope, norms, threshold, method)
    h = ngens_dep(P)
    q = ngens_indep(P)
    c = center(P)
    G = genmat_dep(P)
    GI = genmat_indep(P)
    E = expmat(P)
    idx = indexvector(P)

    # TODO is constructing an array of booleans the most efficient way?
    K = [norms[i] ≤ threshold for i in 1:h]
    Kbar = .!K

    H = [norms[h + i] ≤ threshold for i in 1:q]
    Hbar = .!H

    PZ = SparsePolynomialZonotope(c, G[:, K], GI[:, H], E[:, K], idx)
    Z = reduce_order(overapproximate(PZ, Zonotope), 1, method)

    Ebar = E[:, Kbar]
    N = [!iszero(e) for e in eachrow(Ebar)]

    cz = center(Z)
    Gz = genmat(Z)
    return SparsePolynomialZonotope(cz, G[:, Kbar], hcat(GI[:, Hbar], Gz),
                                    Ebar[N, :], idx[N])
end
