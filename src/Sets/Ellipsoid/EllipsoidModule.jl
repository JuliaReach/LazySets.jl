module EllipsoidModule

using Reexport

using ..LazySets: AbstractCentrallySymmetric, @validate
using LinearAlgebra: dot, I, checksquare, isposdef
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Arrays: inner
using ReachabilityBase.Distribution: reseed!

@reexport import ..API: center, isoperationtype, rand, affine_map, ∈,
                        linear_map, ρ, σ, translate!
@reexport import ..LazySets: ○
@reexport using ..API

export Ellipsoid,
       shape_matrix

include("Ellipsoid.jl")

include("center.jl")
include("in.jl")
include("isoperationtype.jl")
include("rand.jl")
include("affine_map.jl")
include("linear_map.jl")
include("support_function.jl")
include("support_vector.jl")
include("translate.jl")

include("shape_matrix.jl")

function ○(c::VN,
           shape_matrix::MN) where {N<:AbstractFloat,
                                    VN<:AbstractVector{N},
                                    MN<:AbstractMatrix{N}}
    return Ellipsoid(c, shape_matrix)
end

end  # module
