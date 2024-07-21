module HPolyhedronModule

using Reexport, Requires

using ..LazySets: AbstractPolyhedron, LazySet, default_lp_solver,
                  default_polyhedra_backend, iscomplement, is_lp_infeasible,
                  is_lp_optimal, is_lp_unbounded, has_lp_infeasibility_ray,
                  linprog, tosimplehrep, _isempty_polyhedron, _normal_Vector
using ..HalfSpaceModule: HalfSpace
using ..HPolytopeModule: HPolytope
using LinearAlgebra: dot
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Arrays: to_negative_vector
using ReachabilityBase.Basetype: basetype
using ReachabilityBase.Comparison: isapproxzero
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Require: require

@reexport import ..API: constraints_list, dim, isempty, isoperationtype, rand,
                        permute, ρ, σ, translate, convex_hull
@reexport import ..LazySets: is_hyperplanar, normalize,
                             remove_redundant_constraints,
                             remove_redundant_constraints!, tohrep, tovrep,
                             addconstraint!
@reexport import ..Base: convert
@reexport using ..API

export HPolyhedron

include("HPolyhedron.jl")

# convenience union type
const HPoly{N} = Union{HPolytope{N},HPolyhedron{N}}

include("constraints_list.jl")
include("dim.jl")
include("isempty.jl")
include("isoperationtype.jl")
include("rand.jl")
include("permute.jl")
include("support_function.jl")
include("support_vector.jl")
include("translate.jl")
include("convex_hull.jl")

include("is_hyperplanar.jl")
include("normalize.jl")
include("remove_redundant_constraints.jl")
include("tohrep.jl")
include("tovrep.jl")
include("addconstraint.jl")

include("convert.jl")

include("init.jl")

end  # module
