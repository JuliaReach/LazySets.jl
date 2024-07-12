module Line2DModule

using Reexport, Requires

using ..LazySets: AbstractPolyhedron, AbstractLinearMapAlgorithm,
                  _constraints_list_hyperplane, _linear_map_hrep,
                  _non_element_halfspace, _σ_hyperplane_halfspace
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Arrays: nonzero_indices
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Comparison: _isapprox
using ReachabilityBase.Require: require
using LinearAlgebra: dot

@reexport import ..API: an_element, constraints_list, dim, isbounded, isempty,
                        isoperationtype, isuniversal, rand, ∈, project, σ,
                        translate
@reexport import ..LazySets: constrained_dimensions, _linear_map_hrep_helper
@reexport using ..API

export Line2D

include("Line2D.jl")

include("an_element.jl")
include("constrained_dimensions.jl")
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

include("init.jl")

end  # module