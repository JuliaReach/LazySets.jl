module HPolytopeModule

using Reexport, Requires

using ..LazySets: AbstractPolytope, LazySet, AbstractLinearMapAlgorithm,
                  default_polyhedra_backend, vertices_list_1d, _linear_map_hrep,
                  _minkowski_sum_hrep_preprocess, _normal_Vector
using ..HalfSpaceModule: HalfSpace
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Comparison: isapproxzero, _ztol
using ReachabilityBase.Require: require

@reexport import ..API: isbounded, isoperationtype, rand, vertices_list,
                        minkowski_sum
import ..LazySets: _linear_map_hrep_helper, _vertices_list
import Base: convert
@reexport using ..API

export HPolytope

include("HPolytope.jl")

include("isbounded.jl")
include("isoperationtype.jl")
include("rand.jl")
include("vertices_list.jl")
include("linear_map.jl")
include("minkowski_sum.jl")

include("convert.jl")

include("init.jl")

end  # module
