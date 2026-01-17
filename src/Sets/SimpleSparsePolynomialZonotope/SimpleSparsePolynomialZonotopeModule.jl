module SimpleSparsePolynomialZonotopeModule

using Reexport: @reexport

using ..LazySets: AbstractSparsePolynomialZonotope, ngens_dep, nparams,
                  _remove_redundant_generators_polyzono, @validate
using LinearAlgebra: dot, I
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Distribution: reseed!

@reexport import ..API: convex_hull, center, isoperationtype, rand, linear_map,
                        translate, translate!, cartesian_product,
                        linear_combination, minkowski_sum
@reexport import ..LazySets: expmat, genmat, genmat_dep, genmat_indep, ngens,
                             ngens_indep, polynomial_order,
                             remove_redundant_generators
import Base: convert
@reexport using ..API

export SimpleSparsePolynomialZonotope,
       SSPZ,
       quadratic_map

include("SimpleSparsePolynomialZonotope.jl")

# short-hand
const SSPZ = SimpleSparsePolynomialZonotope

include("center.jl")
include("convex_hull.jl")
include("isoperationtype.jl")
include("rand.jl")
include("linear_map.jl")
include("translate.jl")
include("cartesian_product.jl")
include("linear_combination.jl")
include("minkowski_sum.jl")

include("genmat.jl")
include("genmat_dep.jl")
include("genmat_indep.jl")
include("ngens.jl")
include("ngens_indep.jl")
include("polynomial_order.jl")
include("remove_redundant_generators.jl")
include("expmat.jl")
include("quadratic_map.jl")

include("convert.jl")

end  # module
