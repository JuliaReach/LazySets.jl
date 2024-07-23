module TetrahedronModule

using Reexport, Requires

using ..LazySets: AbstractPolytope
using LinearAlgebra: dot, cross
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Comparison: isapproxzero
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Require: require

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

include("init.jl")

end  # module
