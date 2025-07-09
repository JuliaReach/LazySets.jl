module MatrixZonotopeModule

using Reexport

using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Distribution: reseed!
using ..LazySets.SparsePolynomialZonotopeModule
import ..LazySets: linear_map

@reexport import ..API: center, scale, scale!, rand
@reexport import ..LazySets: generators, ngens, indexvector

export MatrixZonotope

include("MatrixZonotope.jl")
include("indexvector.jl")
include("scale.jl")
include("center.jl")
include("rand.jl")

include("generators.jl")
include("ngens.jl")

end #module
