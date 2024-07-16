module VPolygonModule

using Reexport, Requires

using ..LazySets: AbstractPolygon, AbstractHPolygon, HPolygon, halfspace_left,
                  is_right_turn, _area_vlist, _linear_map_vrep
using Random: AbstractRNG, GLOBAL_RNG, shuffle
using ReachabilityBase.Arrays: isabove, rand_pos_neg_zerosum_vector
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Require: require
using LinearAlgebra: dot

@reexport import ..API: an_element, area, constraints_list, isoperationtype,
                        rand, vertices_list, ∈, linear_map, permute, project,
                        σ, translate, translate!
@reexport import ..LazySets: remove_redundant_vertices,
                             remove_redundant_vertices!, tohrep, tovrep
@reexport using ..API

export VPolygon

include("VPolygon.jl")

include("an_element.jl")
include("area.jl")
include("constraints_list.jl")
include("isoperationtype.jl")
include("rand.jl")
include("vertices_list.jl")
include("in.jl")
include("linear_map.jl")
include("permute.jl")
include("project.jl")
include("support_vector.jl")
include("translate.jl")

include("remove_redundant_vertices.jl")
include("tohrep.jl")
include("tovrep.jl")

include("init.jl")

end  # module
