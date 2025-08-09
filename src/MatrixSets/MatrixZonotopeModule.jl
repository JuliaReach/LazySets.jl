module MatrixZonotopeModule

using Reexport

using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Distribution: reseed!
using ..LazySets: Zonotope, genmat

@reexport import ..API: center, scale, scale!, rand, norm, linear_map, minkowski_sum
@reexport import ..LazySets: generators, ngens, order, remove_redundant_generators,
                             reduce_order, AbstractReductionMethod

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

include("generators.jl")
include("indexvector.jl")
include("ngens.jl")
include("order.jl")

end  # module
