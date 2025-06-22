module EmptySetModule

using Reexport, Requires

using ..LazySets: LazySet, ConvexSet, _witness_result_empty
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Commutative: @commutative
using ReachabilityBase.Comparison: _rtol
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Iteration: EmptyIterator
using ReachabilityBase.Require: require

@reexport import ..API: an_element, area, complement, diameter, dim,
                        exponential_map, high, ∈, is_interior_point, isbounded,
                        isboundedtype, isempty, isoperationtype, isuniversal,
                        linear_map, low, norm, permute, project, radius, rand,
                        rectify, reflect, sample, scale, scale!, ρ, σ,
                        translate, translate!, vertices, vertices_list, volume,
                        cartesian_product, convex_hull, difference, distance,
                        intersection, ≈, isdisjoint, ⊆, linear_combination,
                        minkowski_difference, minkowski_sum
@reexport import ..LazySets: chebyshev_center_radius, constrained_dimensions,
                             linear_map_inverse, rationalize, triangulate
import ..LazySets: plot_recipe
import Base: convert, copy
@reexport using ..API

export EmptySet, ∅

include("EmptySet.jl")

include("an_element.jl")
include("area.jl")
include("complement.jl")
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
include("rectify.jl")
include("reflect.jl")
include("sample.jl")
include("rationalize.jl")
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
include("issubset.jl")
include("linear_combination.jl")
include("minkowski_difference.jl")
include("minkowski_sum.jl")

include("chebyshev_center_radius.jl")
include("constrained_dimensions.jl")
include("triangulate.jl")

include("convert.jl")
include("copy.jl")

"""
    plot_recipe(∅::EmptySet{N}, [ε]=zero(N)) where {N}

Convert an empty set to a sequence of points for plotting.
In the special case of an empty set, the sequence is empty.

### Input

- `∅` -- empty set
- `ε` -- (optional, default: `0`) ignored, used for dispatch

### Output

An empty array.
"""
function plot_recipe(::EmptySet{N}, ε=zero(N)) where {N}
    return []
end

include("init.jl")

end  # module
