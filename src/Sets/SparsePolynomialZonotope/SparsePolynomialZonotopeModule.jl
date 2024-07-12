module SparsePolynomialZonotopeModule

using Reexport, Requires

using ..LazySets: AbstractSparsePolynomialZonotope, AbstractReductionMethod,
                  genmat, GIR05, order, remove_zero_columns, Zonotope,
                  _extrema_lowhigh, _remove_redundant_generators_polyzono
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Require: require
import IntervalArithmetic as IA

@reexport import ..API: center, extrema, isoperationtype, rand, linear_map, ρ,
                        translate
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
include("expmat.jl")
include("extrema.jl")
include("genmat_dep.jl")
include("genmat_indep.jl")
include("isoperationtype.jl")
include("polynomial_order.jl")
include("rand.jl")
include("remove_redundant_generators.jl")
include("linear_map.jl")
include("reduce_order.jl")
include("support_function.jl")
include("translate.jl")

include("indexvector.jl")

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
