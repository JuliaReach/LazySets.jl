module StarModule

using Reexport: @reexport

using ..LazySets: AbstractPolyhedron, @validate
using LinearAlgebra: I, dot
using ReachabilityBase.Arrays: At_mul_B, to_matrix

@reexport import ..API: an_element, center, dim, isempty, isoperationtype,
                        affine_map, in, linear_map, ρ, σ
import Base: convert
@reexport using ..API

export Star,
       basis,
       predicate

include("Star.jl")

include("an_element.jl")
include("center.jl")
# include("constraints_list.jl")
include("dim.jl")
# include("isbounded.jl")
include("isempty.jl")
include("isoperationtype.jl")
# include("rand.jl")
# include("vertices_list.jl")
include("affine_map.jl")
include("in.jl")
include("linear_map.jl")
include("support_function.jl")
include("support_vector.jl")

include("basis.jl")
include("predicate.jl")

include("convert.jl")

end  # module
