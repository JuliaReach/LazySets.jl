module MatrixZonotopeModule

using Reexport

using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Commutative
using ..LazySets.SparsePolynomialZonotopeModule
using ..LazySets: Zonotope

@reexport import ..API: center, scale, scale!, rand, norm, linear_map
@reexport import ..LazySets: generators, ngens, indexvector

export AbstractMatrixZonotope, MatrixZonotope, MatrixZonotopeProduct,
       factors, nfactors, remove_redundant_factors

include("AbstractMatrixZonotope.jl")
include("MatrixZonotope.jl")
include("MatrixZonotopeProduct.jl")
include("indexvector.jl")
include("norm.jl")
include("linear_map.jl")
include("scale.jl")
include("center.jl")
include("rand.jl")
include("scale.jl")
include("norm.jl")

include("generators.jl")
include("indexvector.jl")
include("ngens.jl")

end #module
