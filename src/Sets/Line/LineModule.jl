module LineModule

using Reexport, Requires

using ..LazySets: AbstractPolyhedron, @validate, @validate_commutative
import LinearAlgebra
using LinearAlgebra: dot, nullspace
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Arrays: ismultiple
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Comparison: _isapprox, isapproxzero
using ReachabilityBase.Require

@reexport import ..API: an_element, constraints_list, dim, isbounded, isempty,
                        isoperationtype, isuniversal, project, rand, distance,
                        ∈, linear_map, ρ, σ, translate!
@reexport import ..LazySets: normalize
@reexport import LinearAlgebra: normalize!
@reexport using ..API

export Line,
       direction

include("Line.jl")

include("an_element.jl")
include("constraints_list.jl")
include("dim.jl")
include("isbounded.jl")
include("isempty.jl")
include("isoperationtype.jl")
include("isuniversal.jl")
include("normalize.jl")
include("project.jl")
include("rand.jl")
include("distance.jl")
include("in.jl")
include("linear_map.jl")
include("support_function.jl")
include("support_vector.jl")
include("translate.jl")

include("direction.jl")

include("init.jl")

end  # module
