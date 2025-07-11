module MatrixZonotopeModule

using Reexport

using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Distribution: reseed!
using ..LazySets.SparsePolynomialZonotopeModule
using ..LazySets: Zonotope

@reexport import ..API: center, scale, scale!, rand, norm
@reexport import ..LazySets: generators, ngens, indexvector

export MatrixZonotope

include("MatrixZonotope.jl")
include("indexvector.jl")
include("norm.jl")
include("scale.jl")
include("center.jl")
include("rand.jl")

include("generators.jl")
include("ngens.jl")

end #module
