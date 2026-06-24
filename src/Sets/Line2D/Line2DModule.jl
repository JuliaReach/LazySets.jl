module Line2DModule

using Reexport: @reexport
using ..LazySets: AbstractPolyhedron, _intersection_line2d,
                  _witness_result_empty, @validate
using LinearAlgebra: dot
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Arrays: nonzero_indices, right_turn
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Comparison: isapproxzero, _isapprox

@reexport import ..API: an_element, dim, isbounded, isempty, isoperationtype,
                        rand, in, project, scale, scale!, translate,
                        intersection, isdisjoint
@reexport import ..LazySets: constrained_dimensions
@reexport using ..API

export Line2D

include("Line2D.jl")

include("an_element.jl")
# include("constraints_list.jl")
include("dim.jl")
include("isbounded.jl")
include("isempty.jl")
include("isoperationtype.jl")
# include("isuniversal.jl")
include("rand.jl")
include("in.jl")
# include("linear_map.jl")
include("project.jl")
include("scale.jl")
# include("support_vector.jl")
include("translate.jl")
include("intersection.jl")
include("isdisjoint.jl")

include("constrained_dimensions.jl")

end  # module
