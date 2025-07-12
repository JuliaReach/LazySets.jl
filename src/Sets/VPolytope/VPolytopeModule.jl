module VPolytopeModule

using Reexport, Requires

using ..LazySets: AbstractPolytope, LazySet, LinearMapVRep, default_lp_solver,
                  default_lp_solver_polyhedra, default_polyhedra_backend,
                  is_lp_infeasible, is_lp_optimal, linprog, _extrema_vlist,
                  _high_vlist, _low_vlist, _minkowski_sum_vrep_nd,
                  _removevredundancy!, _scale_copy_inplace, @validate
using LinearAlgebra: dot
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Arrays: projection_matrix
using ReachabilityBase.Comparison: _ztol
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Require: require

@reexport import ..API: constraints_list, dim, extrema, high, isoperationtype,
                        low, rand, reflect, vertices_list, ∈, linear_map,
                        permute, project, scale, scale!, ρ, σ, translate,
                        translate!, cartesian_product, convex_hull,
                        minkowski_sum
@reexport import ..LazySets: remove_redundant_vertices, tohrep, tovrep
import ..LazySets: _linear_map_vrep
import Base: convert
@reexport using ..API

export VPolytope

include("VPolytope.jl")

include("constraints_list.jl")
include("dim.jl")
include("extrema.jl")
include("high.jl")
include("isoperationtype.jl")
include("low.jl")
include("rand.jl")
include("reflect.jl")
include("vertices_list.jl")
include("in.jl")
include("linear_map.jl")
include("permute.jl")
include("project.jl")
include("scale.jl")
include("support_function.jl")
include("support_vector.jl")
include("translate.jl")
include("cartesian_product.jl")
include("convex_hull.jl")
include("minkowski_sum.jl")

include("polyhedron.jl")
include("remove_redundant_vertices.jl")
include("tohrep.jl")
include("tovrep.jl")

include("convert.jl")

include("init.jl")

end  # module
