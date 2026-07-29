module HyperplaneModule

using Reexport: @reexport

using ..LazySets: AbstractPolyhedron, @validate
using LinearAlgebra: dot
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Arrays: nonzero_indices
using ReachabilityBase.Commutative: @commutative
using ReachabilityBase.Comparison: _isapprox
using ReachabilityBase.Distribution: reseed!

@reexport import ..API: an_element, constraints_list, dim, isbounded, isempty,
                        isoperationtype, rand, reflect, in, project, ρ, σ,
                        translate, isdisjoint
@reexport import ..LazySets: constrained_dimensions, ishyperplanar
@reexport using ..API

export Hyperplane

include("Hyperplane.jl")

include("an_element.jl")
include("constraints_list.jl")
include("dim.jl")
include("isbounded.jl")
include("isempty.jl")
include("isoperationtype.jl")
# include("isuniversal.jl")
include("rand.jl")
include("reflect.jl")
# include("distance.jl")
include("in.jl")
include("linear_map.jl")
include("project.jl")
include("support_function.jl")
include("support_vector.jl")
include("translate.jl")
include("isdisjoint.jl")

include("constrained_dimensions.jl")
include("ishyperplanar.jl")
# include("normalize.jl")

# include("convert.jl")

end  # module
