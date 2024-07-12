module SimpleSparsePolynomialZonotopeModule

using Reexport

using ..LazySets: AbstractPolynomialZonotope, nparams,
                  _remove_redundant_generators_polyzono
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Comparison: isapproxzero
using LinearAlgebra: dot

@reexport import ..API: convex_hull, center, isoperationtype, rand, linear_map
@reexport import ..LazySets: expmat, genmat, ngens, polynomial_order,
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
include("isoperationtype.jl")
include("ngens.jl")
include("ngens_dep.jl")
include("ngens_indep.jl")
include("polynomial_order.jl")
include("remove_redundant_generators.jl")
include("expmat.jl")
include("quadratic_map.jl")
include("rand.jl")
include("linear_map.jl")

end  # module
