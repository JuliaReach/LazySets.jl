module UniverseModule

using Reexport, Requires

using ..LazySets: LazySet, AbstractPolyhedron, default_polyhedra_backend,
                  _witness_result_empty
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Commutative: @commutative
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Iteration: EmptyIterator
using ReachabilityBase.Require: require

@reexport import ..API: an_element, area, complement, constraints,
                        constraints_list, diameter, dim, isbounded,
                        isboundedtype, isempty, isoperationtype, isuniversal,
                        norm, radius, rand, reflect, volume, ∈, permute,
                        project, scale, scale!, ρ, σ, translate, translate!,
                        cartesian_product, convex_hull, difference, distance,
                        intersection, isdisjoint, ⊆, linear_combination,
                        minkowski_difference, minkowski_sum
@reexport import ..LazySets: constrained_dimensions, linear_map_inverse,
                             rationalize, tosimplehrep, triangulate
import Base: copy
@reexport using ..API

export Universe

include("Universe.jl")

include("an_element.jl")
include("area.jl")
include("complement.jl")
include("constraints.jl")
include("constraints_list.jl")
include("copy.jl")
include("diameter.jl")
include("dim.jl")
include("isbounded.jl")
include("isboundedtype.jl")
include("isdisjoint.jl")
include("isempty.jl")
include("isoperationtype.jl")
include("isuniversal.jl")
include("norm.jl")
include("radius.jl")
include("rand.jl")
include("rationalize.jl")
include("reflect.jl")
include("volume.jl")
include("in.jl")
include("permute.jl")
include("project.jl")
include("scale.jl")
include("support_function.jl")
include("support_vector.jl")
include("translate.jl")
include("cartesian_product.jl")
include("convex_hull.jl")
include("difference.jl")
include("distance.jl")
include("intersection.jl")
include("issubset.jl")
include("linear_combination.jl")
include("minkowski_difference.jl")
include("minkowski_sum.jl")

include("constrained_dimensions.jl")
include("linear_map_inverse.jl")
include("tosimplehrep.jl")
include("triangulate.jl")

include("init.jl")

end  # module
