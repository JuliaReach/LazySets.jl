module SparsePolynomialZonotopeModule

using Reexport, Requires

using ..LazySets: AbstractSparsePolynomialZonotope, AbstractReductionMethod,
                  genmat, GIR05, order, _remove_redundant_generators_polyzono,
                  MatrixZonotope, MatrixZonotopeProduct, ngens, generators,
                  factors, @validate
import IntervalArithmetic as IA
using LinearAlgebra: I
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Arrays: remove_zero_columns
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Require: require

@reexport import ..API: center, isoperationtype, rand, scale, translate,
                        translate!, exact_sum, linear_map
@reexport import ..LazySets: expmat, genmat_dep, genmat_indep, indexvector,
                             ngens_dep, ngens_indep, nparams, polynomial_order,
                             reduce_order, remove_redundant_generators
@reexport using ..API

export SparsePolynomialZonotope

include("SparsePolynomialZonotope.jl")

# short-hand
const SPZ = SparsePolynomialZonotope

include("center.jl")
include("isoperationtype.jl")
include("rand.jl")
include("scale.jl")
include("translate.jl")
include("exact_sum.jl")
include("linear_map.jl")

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
uniqueID(n::Int) = collect(1:n)

end  # module