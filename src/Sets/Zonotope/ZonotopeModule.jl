module ZonotopeModule

using Reexport, Requires

using ..LazySets: AbstractZonotope, generators_fallback,
                  _vertices_list_zonotope_iterative
using LinearAlgebra: mul!
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Arrays: ismultiple, remove_zero_columns, to_matrix,
                               vector_type
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Require: require

@reexport import ..API: center, high, isoperationtype, low, rand, vertices_list,
                        permute, scale!, translate!
@reexport import ..LazySets: generators, genmat, ngens, reduce_order,
                             remove_redundant_generators, togrep
import Base: convert
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
include("vertices_list.jl")
include("linear_map.jl")
include("permute.jl")
include("scale.jl")
include("translate.jl")

include("generators.jl")
include("genmat.jl")
include("ngens.jl")
include("reduce_order.jl")
include("remove_redundant_generators.jl")
include("togrep.jl")
include("remove_zero_generators.jl")
include("split.jl")

include("convert.jl")

include("init.jl")

end  # module
