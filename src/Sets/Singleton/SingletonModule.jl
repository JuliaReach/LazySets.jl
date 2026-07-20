module SingletonModule

using Reexport: @reexport

using ..LazySets: AbstractSingleton, @validate
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Distribution: reseed!

@reexport import ..API: center, isoperationtype, rand, rectify, linear_map,
                        permute, project, scale, scale!, translate!
@reexport import ..LazySets: singleton_list
import Base: copy
@reexport using ..API

export Singleton

include("Singleton.jl")

include("center.jl")
include("isoperationtype.jl")
include("rand.jl")
include("rectify.jl")
include("linear_map.jl")
include("permute.jl")
include("project.jl")
include("scale.jl")
include("translate.jl")

include("singleton_list.jl")

include("copy.jl")

end  # module
