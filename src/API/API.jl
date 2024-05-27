"""
    API

This module contains an API (application programming interface) for set
libraries. The module only defines and documents the general functions and does
not provide implementations.
"""
module API

import Base: eltype, extrema, isdisjoint, isempty, ∈, ≈, ==, ⊆
import Random: rand
import LinearAlgebra: norm
import SparseArrays: permute
import ReachabilityBase.Arrays: distance, rectify

export
# unary set operations
      an_element, area, center, complement, concretize, constraints_list,
      constraints, convex_hull, diameter, dim, high, is_polyhedral, isbounded,
      isboundedtype, isconvextype, isempty, isoperation, isoperationtype,
      isuniversal, low, norm, radius, rectify, reflect, sample,
      support_function, ρ, support_vector, σ, surface, vertices_list, vertices,
      volume,
# mixed set operations (typically with vectors or matrices)
      affine_map, exponential_map, is_interior_point, linear_map, project,
      scale!, scale, support_function, ρ, support_vector, σ, translate!,
      translate,
# binary set operations
      cartesian_product, difference, distance, exact_sum, ⊞, intersection,
      is_intersection_empty, isequivalent, ⊂, linear_combination,
      minkowski_difference, pontryagin_difference, minkowski_sum

include("LazySet.jl")

include("Unary/an_element.jl")
include("Unary/area.jl")
include("Unary/center.jl")
include("Unary/complement.jl")
include("Unary/concretize.jl")
include("Unary/constraints_list.jl")
include("Unary/constraints.jl")
include("Unary/convex_hull.jl")
include("Unary/diameter.jl")
include("Unary/dim.jl")
include("Unary/eltype.jl")
include("Unary/extrema.jl")
include("Unary/high.jl")
include("Unary/is_polyhedral.jl")
include("Unary/isbounded.jl")
include("Unary/isboundedtype.jl")
include("Unary/isconvextype.jl")
include("Unary/isempty.jl")
include("Unary/isoperation.jl")
include("Unary/isoperationtype.jl")
include("Unary/isuniversal.jl")
include("Unary/low.jl")
include("Unary/norm.jl")
include("Unary/radius.jl")
include("Unary/rand.jl")
include("Unary/rectify.jl")
include("Unary/reflect.jl")
include("Unary/surface.jl")
include("Unary/vertices_list.jl")
include("Unary/vertices.jl")
include("Unary/volume.jl")

include("Mixed/affine_map.jl")
include("Mixed/exponential_map.jl")
include("Mixed/in.jl")
include("Mixed/is_interior_point.jl")
include("Mixed/linear_map.jl")
include("Mixed/permute.jl")
include("Mixed/project.jl")
include("Mixed/sample.jl")
include("Mixed/scale!.jl")
include("Mixed/scale.jl")
include("Mixed/support_function.jl")
include("Mixed/support_vector.jl")
include("Mixed/translate!.jl")
include("Mixed/translate.jl")

include("Binary/cartesian_product.jl")
include("Binary/difference.jl")
include("Binary/distance.jl")
include("Binary/exact_sum.jl")
include("Binary/intersection.jl")
include("Binary/isapprox.jl")
include("Binary/isdisjoint.jl")
include("Binary/isequal.jl")
include("Binary/isequivalent.jl")
include("Binary/isstrictsubset.jl")
include("Binary/issubset.jl")
include("Binary/linear_combination.jl")
include("Binary/minkowski_difference.jl")
include("Binary/minkowski_sum.jl")

end  # module
