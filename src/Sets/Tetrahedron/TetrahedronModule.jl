module TetrahedronModule

using Reexport

using ..LazySets: AbstractPolytope, VPolytope
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Comparison: isapproxzero
using LinearAlgebra: dot, cross

@reexport import ..API: constraints_list, dim, isoperationtype, rand, ∈, σ
@reexport using ..API

export Tetrahedron

include("Tetrahedron.jl")

include("constraints_list.jl")
include("dim.jl")
include("isoperationtype.jl")
include("rand.jl")
include("in.jl")
include("support_vector.jl")

end  # module
