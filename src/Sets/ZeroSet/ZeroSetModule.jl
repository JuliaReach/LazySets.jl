module ZeroSetModule

using Reexport: @reexport

using ..LazySets: AbstractSingleton, @validate
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Distribution: reseed!

@reexport import ..API: dim, isoperationtype, rand, rectify, reflect, in,
                        linear_map, scale, scale!, ρ, minkowski_sum
@reexport import ..LazySets: element
import Base: copy
@reexport using ..API

export ZeroSet

include("ZeroSet.jl")

include("copy.jl")
include("dim.jl")
include("isoperationtype.jl")
include("rand.jl")
include("rectify.jl")
include("reflect.jl")
include("in.jl")
include("linear_map.jl")
include("scale.jl")
include("support_function.jl")
# include("translate.jl")
include("minkowski_sum.jl")

include("element.jl")

end  # module
