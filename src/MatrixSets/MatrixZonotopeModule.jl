module MatrixZonotopeModule

using Reexport: @reexport
using Requires: @require

using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Distribution: reseed!
using ..LazySets: Zonotope, genmat, AbstractReductionMethod, dim, GIR05

@reexport import ..API: center, scale, scale!, rand, norm, linear_map, minkowski_sum
@reexport import ..LazySets: generators, ngens, order, remove_redundant_generators,
                             reduce_order

export AbstractMatrixZonotope, MatrixZonotope, MatrixZonotopeProduct,
       MatrixZonotopeExp, indexvector, factors, nfactors

include("AbstractMatrixZonotope.jl")
include("MatrixZonotope.jl")
include("MatrixZonotopeProduct.jl")
include("MatrixZonotopeExp.jl")
include("center.jl")
include("norm.jl")
include("rand.jl")
include("linear_map.jl")
include("scale.jl")
include("minkowski_sum.jl")
include("reduce_order.jl")
include("remove_redundant_generators.jl")
include("reshape.jl")
include("convert.jl")

include("generators.jl")
include("indexvector.jl")
include("ngens.jl")
include("order.jl")

include("init.jl")

end  # module
