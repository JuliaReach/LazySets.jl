module SimpleSparsePolynomialZonotopeModule

using Reexport

using ..LazySets: AbstractSparsePolynomialZonotope, ngens_dep, nparams,
                  _remove_redundant_generators_polyzono
using LinearAlgebra: dot
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Comparison: isapproxzero

@reexport import ..API: convex_hull, center, isoperationtype, rand, linear_map
@reexport import ..LazySets: expmat, genmat, genmat_dep, genmat_indep, ngens,
                             ngens_indep, polynomial_order,
                             remove_redundant_generators
@reexport using ..API

export SimpleSparsePolynomialZonotope,
       SSPZ,
       quadratic_map

include("SimpleSparsePolynomialZonotope.jl")

# short-hand
const SSPZ = SimpleSparsePolynomialZonotope

include("center.jl")
include("convex_hull.jl")
include("genmat.jl")
include("genmat_dep.jl")
include("genmat_indep.jl")
include("isoperationtype.jl")
include("ngens.jl")
include("ngens_indep.jl")
include("polynomial_order.jl")
include("remove_redundant_generators.jl")
include("expmat.jl")
include("quadratic_map.jl")
include("rand.jl")
include("linear_map.jl")

end  # module
