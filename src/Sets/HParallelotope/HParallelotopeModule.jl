module HParallelotopeModule

using Reexport, Requires

using ..LazySets: AbstractZonotope, HalfSpace, generators_fallback, order,
                  _constraints_list_zonotope
using LinearAlgebra: checksquare, det
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Arrays: to_negative_vector
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Require: require

@reexport import ..API: center, constraints_list, dim, isoperationtype, rand,
                        volume
@reexport import ..LazySets: generators, genmat
import Base: convert
@reexport using ..API

export HParallelotope,
       directions,
       base_vertex,
       extremal_vertices,
       offset

include("HParallelotope.jl")

include("base_vertex.jl")
include("center.jl")
include("constraints_list.jl")
include("dim.jl")
include("directions.jl")
include("extremal_vertices.jl")
include("generators.jl")
include("genmat.jl")
include("isoperationtype.jl")
include("offset.jl")
include("rand.jl")
include("volume.jl")

include("convert.jl")

include("init.jl")

end  # module
