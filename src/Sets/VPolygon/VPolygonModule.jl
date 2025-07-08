module VPolygonModule

using Reexport, Requires

using ..LazySets: AbstractPolygon, LazySet, AbstractHPolygon, halfspace_left,
                  is_right_turn, _area_vlist_2D, _extrema_vlist, _high_vlist,
                  _infeasible_constraints_list, _intersection_vrep_2d,
                  _linear_map_vrep, _low_vlist, _minkowski_sum_vrep_2d,
                  _to_colVector, @validate
using ..HPolygonModule: HPolygon
using LinearAlgebra: dot
using Random: AbstractRNG, GLOBAL_RNG, shuffle
using ReachabilityBase.Arrays: isabove, rand_pos_neg_zerosum_vector
using ReachabilityBase.Comparison: _isapprox
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Require: require

@reexport import ..API: an_element, area, constraints_list, extrema, high,
                        isoperationtype, low, rand, vertices_list, ∈,
                        linear_map, permute, project, σ, translate, translate!,
                        convex_hull, intersection, minkowski_sum
@reexport import ..LazySets: remove_redundant_vertices,
                             remove_redundant_vertices!, tohrep, tovrep
import Base: convert
@reexport using ..API

export VPolygon

include("VPolygon.jl")

include("an_element.jl")
include("area.jl")
include("constraints_list.jl")
include("extrema.jl")
include("high.jl")
include("isoperationtype.jl")
include("low.jl")
include("rand.jl")
include("vertices_list.jl")
include("in.jl")
include("linear_map.jl")
include("permute.jl")
include("project.jl")
include("support_vector.jl")
include("translate.jl")
include("convex_hull.jl")
include("intersection.jl")
include("minkowski_sum.jl")

include("remove_redundant_vertices.jl")
include("tohrep.jl")
include("tovrep.jl")

include("convert.jl")

include("init.jl")

end  # module
