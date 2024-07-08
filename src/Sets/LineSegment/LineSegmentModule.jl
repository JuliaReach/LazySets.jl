module LineSegmentModule

using Reexport

using ..LazySets: AbstractZonotope, right_turn
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Comparison: _isapprox, isapproxzero, _leq
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Iteration: EmptyIterator, SingletonIterator
using LinearAlgebra: dot

@reexport import ..API: an_element, center, constraints_list, dim,
                        isoperationtype, rand, vertices_list, ∈, scale!, ρ, σ,
                        translate
@reexport import ..LazySets: generators, genmat, ngens, halfspace_left,
                             halfspace_right
@reexport using ..API

export LineSegment

include("LineSegment.jl")

include("an_element.jl")
include("center.jl")
include("constraints_list.jl")
include("dim.jl")
include("generators.jl")
include("genmat.jl")
include("isoperationtype.jl")
include("ngens.jl")
include("rand.jl")
include("vertices_list.jl")
include("in.jl")
include("scale.jl")
include("support_function.jl")
include("support_vector.jl")
include("translate.jl")

include("halfspace_left.jl")
include("halfspace_right.jl")

end  # module
