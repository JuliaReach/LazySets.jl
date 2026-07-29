module EmptySetModule

using Reexport: @reexport

using ..LazySets: LazySet, ConvexSet, _witness_result, _witness_result_empty,
                  @validate, @validate_commutative
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Comparison: _rtol
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Iteration: EmptyIterator

@reexport import ..API: an_element, area, diameter, dim, exponential_map, high,
                        in, is_interior_point, isbounded, isboundedtype,
                        isempty, isoperationtype, isuniversal, linear_map, low,
                        norm, permute, project, radius, rand, rectify, reflect,
                        sample, scale, scale!, ρ, σ, translate, translate!,
                        vertices, vertices_list, volume, cartesian_product,
                        convex_hull, difference, distance, intersection,
                        isapprox, isdisjoint, isequivalent, ⊂, issubset,
                        linear_combination, minkowski_difference, minkowski_sum
@reexport import ..LazySets: chebyshev_center_radius, constrained_dimensions,
                             linear_map_inverse, rationalize, triangulate
import Base: convert, copy
@reexport using ..API

export EmptySet, ∅

include("EmptySet.jl")

include("an_element.jl")
include("area.jl")
# include("complement.jl")
include("diameter.jl")
include("dim.jl")
include("high.jl")
include("isbounded.jl")
include("isboundedtype.jl")
include("isempty.jl")
include("isoperationtype.jl")
include("isuniversal.jl")
include("low.jl")
include("norm.jl")
include("radius.jl")
include("rand.jl")
include("rationalize.jl")
include("rectify.jl")
include("reflect.jl")
include("sample.jl")
include("vertices_list.jl")
include("vertices.jl")
include("volume.jl")
include("exponential_map.jl")
include("in.jl")
include("is_interior_point.jl")
include("linear_map.jl")
include("linear_map_inverse.jl")
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
include("isapprox.jl")
include("isdisjoint.jl")
include("isequivalent.jl")
include("isstrictsubset.jl")
include("issubset.jl")
include("linear_combination.jl")
include("minkowski_difference.jl")
include("minkowski_sum.jl")

include("chebyshev_center_radius.jl")
include("constrained_dimensions.jl")
# include("polyhedron.jl")
include("triangulate.jl")

include("convert.jl")
include("copy.jl")

end  # module
