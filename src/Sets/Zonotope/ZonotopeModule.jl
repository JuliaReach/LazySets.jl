module ZonotopeModule

using Reexport: @reexport
using Requires: @require

using ..LazySets: AbstractZonotope, generators_fallback, @validate
using LinearAlgebra: mul!
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Arrays: remove_zero_columns, to_matrix
using ReachabilityBase.Distribution: reseed!

@reexport import ..API: center, high, isoperationtype, low, rand, affine_map,
                        permute, scale, scale!, translate, translate!
@reexport import ..LazySets: generators, genmat, ngens, rationalize,
                             reduce_order, remove_redundant_generators
import Base: convert, copy
@reexport using ..API

export Zonotope,
       remove_zero_generators,
       linear_map!,
       split!

include("Zonotope.jl")

include("center.jl")
include("high.jl")
include("isoperationtype.jl")
include("low.jl")
include("rand.jl")
include("rationalize.jl")
include("affine_map.jl")
include("linear_map.jl")
include("permute.jl")
include("scale.jl")
include("translate.jl")

include("generators.jl")
include("genmat.jl")
include("ngens.jl")
include("reduce_order.jl")
include("remove_redundant_generators.jl")
include("remove_zero_generators.jl")
include("split.jl")

include("convert.jl")
include("copy.jl")

include("init.jl")

end  # module
