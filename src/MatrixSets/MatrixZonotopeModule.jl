module MatrixZonotopeModule

using Reexport

using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Distribution: reseed!
using ..LazySets.SparsePolynomialZonotopeModule
using ..LazySets: Zonotope

@reexport import ..API: center, scale, scale!, rand, norm
@reexport import ..LazySets: generators, ngens
@reexport import ..LazySets.SparsePolynomialZonotopeModule: indexvector

export MatrixZonotope

include("MatrixZonotope.jl")

include("center.jl")
include("rand.jl")
include("scale.jl")
include("norm.jl")

include("generators.jl")
include("indexvector.jl")
include("ngens.jl")

end #module
