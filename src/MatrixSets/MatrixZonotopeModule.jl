module MatrixZonotopeModule

using Reexport

using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Distribution: reseed!
using ..LazySets.SparsePolynomialZonotopeModule

@reexport import ..API: center, scale, scale!, rand, norm
@reexport import ..LazySets: generators, ngens, indexvector

export MatrixZonotope, _matrixzonotope_norm

include("MatrixZonotope.jl")
include("indexvector.jl")
include("norm.jl")
include("scale.jl")
include("center.jl")
include("rand.jl")

include("generators.jl")
include("ngens.jl")

end #module
