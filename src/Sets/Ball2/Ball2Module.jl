module Ball2Module

using Reexport: @reexport
using Requires: @require

using ..LazySets: AbstractBallp, _witness_result_empty, @validate
using LinearAlgebra: dot, axpby!
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Comparison: _leq, isapproxzero
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Require: require

@reexport import ..API: area, center, isoperationtype, rand, reflect, sample,
                        volume, in, project, scale, ρ, σ, translate!,
                        isdisjoint, issubset
@reexport import ..LazySets: ○, chebyshev_center_radius, norm_ball, radius_ball
@reexport using ..API

export Ball2

include("Ball2.jl")

include("area.jl")
include("center.jl")
include("in.jl")
include("isoperationtype.jl")
include("project.jl")
include("rand.jl")
include("reflect.jl")
include("sample.jl")
include("scale.jl")
include("support_function.jl")
include("support_vector.jl")
include("translate.jl")
include("volume.jl")
include("isdisjoint.jl")
include("issubset.jl")

include("chebyshev_center_radius.jl")
include("norm_ball.jl")
include("radius_ball.jl")

function ○(c::VN, r::N) where {N<:AbstractFloat,VN<:AbstractVector{N}}
    return Ball2(c, r)
end

include("init.jl")

end  # module
