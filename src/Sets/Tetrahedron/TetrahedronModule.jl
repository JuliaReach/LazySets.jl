module TetrahedronModule

using Reexport: @reexport

using ..LazySets: AbstractPolytope, @validate
using LinearAlgebra: dot, cross
using ReachabilityBase.Comparison: isapproxzero

@reexport import ..API: dim, isoperationtype, in
@reexport using ..API

export Tetrahedron

include("Tetrahedron.jl")

# include("constraints_list.jl")
include("dim.jl")
include("isoperationtype.jl")
# include("rand.jl")
include("in.jl")
# include("support_vector.jl")

end  # module
