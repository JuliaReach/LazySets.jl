module SparsePolynomialZonotopeModule

using Reexport: @reexport

using ..LazySets: AbstractSparsePolynomialZonotope, AbstractReductionMethod,
                  GIR05, order, _remove_redundant_generators_polyzono,
                  MatrixZonotope, MatrixZonotopeProduct, ngens, generators,
                  factors, @validate
import IntervalArithmetic as IA
import LazySets: _indexvector
using LinearAlgebra: I

@reexport import ..API: center, isoperationtype, scale, scale!, translate,
                        translate!, exact_sum, linear_map
@reexport import ..LazySets: expmat, genmat_dep, genmat_indep, indexvector,
                             ngens_dep, ngens_indep, nparams, polynomial_order,
                             reduce_order, remove_redundant_generators
@reexport using ..API

export SparsePolynomialZonotope

include("SparsePolynomialZonotope.jl")

include("center.jl")
include("isoperationtype.jl")
# include("rand.jl")
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

"""
    uniqueID(n::Int)

Return a collection of n unique identifiers (integers 1, …, n).

### Input

- `n` -- number of variables

### Output

`1:n`.
"""
uniqueID(n::Int) = collect(1:n)
_indexvector(P::SparsePolynomialZonotope) = indexvector(P)

end  # module
