module ZeroSetModule

using Reexport

using ..LazySets: AbstractSingleton, Singleton
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Distribution: reseed!

@reexport import ..API: dim, isoperationtype, rand, rectify, reflect, ∈,
                        linear_map, scale, scale!, ρ, translate, minkowski_sum
@reexport import ..LazySets: element
@reexport using ..API

export ZeroSet

include("ZeroSet.jl")

include("dim.jl")
include("isoperationtype.jl")
include("rand.jl")
include("rectify.jl")
include("reflect.jl")
include("in.jl")
include("linear_map.jl")
include("scale.jl")
include("support_function.jl")
include("translate.jl")

include("element.jl")

include("minkowski_sum.jl")

end  # module
