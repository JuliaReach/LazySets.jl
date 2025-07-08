module MatrixZonotope

using Reexport, Requires

using ..LazySets.SparsePolynomialZonotopeModule
import ..LazySets: linear_map

@reexport import ..API: center, scale, scale!
@reexport import ..LazySets: generators, ngens

include("MatrixZonotope.jl")
include("center.jl")
include("generators.jl")
include("ngens.jl")
include("scale.jl")

end
