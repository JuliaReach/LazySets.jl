module MatrixZonotopeModule

using Reexport

using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Distribution: reseed!
using ..LazySets: Zonotope

@reexport import ..API: center, scale, scale!, rand, norm, linear_map
@reexport import ..LazySets: generators, ngens

export AbstractMatrixZonotope, MatrixZonotope, MatrixZonotopeProduct,
       MatrixZonotopeExp, indexvector, factors, nfactors

include("AbstractMatrixZonotope.jl")
include("MatrixZonotope.jl")
include("MatrixZonotopeProduct.jl")
include("MatrixZonotopeExp.jl")
include("linear_map.jl")
include("scale.jl")
include("center.jl")
include("rand.jl")
include("norm.jl")

include("generators.jl")
include("indexvector.jl")
include("ngens.jl")

end #module