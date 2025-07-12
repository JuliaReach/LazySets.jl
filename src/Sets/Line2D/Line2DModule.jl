module Line2DModule

using Reexport, Requires

using ..LazySets: AbstractPolyhedron, AbstractLinearMapAlgorithm,
                  _constraints_list_hyperplane, _intersection_line2d,
                  _linear_map_hrep, _non_element_halfspace,
                  _σ_hyperplane_halfspace, _witness_result_empty, @validate
using LinearAlgebra: dot
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Arrays: nonzero_indices, right_turn
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Comparison: isapproxzero, _isapprox
using ReachabilityBase.Require: require

@reexport import ..API: an_element, constraints_list, dim, isbounded, isempty,
                        isoperationtype, isuniversal, rand, ∈, project, σ,
                        translate, intersection, isdisjoint
@reexport import ..LazySets: constrained_dimensions
import ..LazySets: _linear_map_hrep_helper
@reexport using ..API

export Line2D

include("Line2D.jl")

include("an_element.jl")
include("constraints_list.jl")
include("dim.jl")
include("isbounded.jl")
include("isempty.jl")
include("isoperationtype.jl")
include("isuniversal.jl")
include("rand.jl")
include("in.jl")
include("linear_map.jl")
include("project.jl")
include("support_vector.jl")
include("translate.jl")
include("intersection.jl")
include("isdisjoint.jl")

include("constrained_dimensions.jl")

include("init.jl")

end  # module
