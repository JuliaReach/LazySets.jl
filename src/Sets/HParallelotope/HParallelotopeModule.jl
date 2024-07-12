module HParallelotopeModule

using Reexport

using ..LazySets: AbstractZonotope, HalfSpace, HPolyhedron, generators_fallback
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Arrays: to_negative_vector
using ReachabilityBase.Distribution: reseed!
using LinearAlgebra: checksquare, det

@reexport import ..API: center, constraints_list, dim, isoperationtype, rand,
                        volume
@reexport import ..LazySets: generators, genmat
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

end  # module
