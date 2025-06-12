module SparsePolynomialZonotopeModule

using Reexport, Requires

using ..LazySets: AbstractSparsePolynomialZonotope, AbstractReductionMethod,
                  genmat, GIR05, order, _extrema_lowhigh,
                  _remove_redundant_generators_polyzono
import IntervalArithmetic as IA
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Arrays: remove_zero_columns
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Require: require

@reexport import ..API: center, extrema, isoperationtype, rand, linear_map,
                        scale, ρ, translate, exact_sum
@reexport import ..LazySets: expmat, genmat_dep, genmat_indep, ngens_dep,
                             ngens_indep, nparams, polynomial_order,
                             reduce_order, remove_redundant_generators
@reexport using ..API

export SparsePolynomialZonotope,
       indexvector

include("SparsePolynomialZonotope.jl")

# short-hand
const SPZ = SparsePolynomialZonotope

include("center.jl")
include("extrema.jl")
include("isoperationtype.jl")
include("rand.jl")
include("linear_map.jl")
include("scale.jl")
include("support_function.jl")
include("translate.jl")
include("exact_sum.jl")

include("expmat.jl")
include("genmat_dep.jl")
include("genmat_indep.jl")
include("indexvector.jl")
include("merge_id.jl")
include("polynomial_order.jl")
include("remove_redundant_generators.jl")
include("reduce_order.jl")

include("init.jl")

"""
    uniqueID(n::Int)

Return a collection of n unique identifiers (integers 1, …, n).

### Input

- `n` -- number of variables

### Output

`1:n`.
"""
uniqueID(n::Int) = 1:n

end  # module
