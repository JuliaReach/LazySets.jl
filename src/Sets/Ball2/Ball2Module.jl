module Ball2Module

using Reexport, Requires

using ..LazySets: AbstractBallp, _witness_result_empty
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Comparison: _leq, isapproxzero
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Require
using LinearAlgebra: dot, axpby!

@reexport import ..API: area, center, isoperationtype, rand, reflect, sample,
                        volume, ∈, project, scale, ρ, σ, translate!,
                        isdisjoint, ⊆
@reexport import ..LazySets: ball_norm, chebyshev_center_radius, ○, radius_ball
@reexport using ..API

export Ball2

include("Ball2.jl")

include("area.jl")
include("ball_norm.jl")
include("center.jl")
include("chebyshev_center_radius.jl")
include("in.jl")
include("isoperationtype.jl")
include("project.jl")
include("radius_ball.jl")
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

function ○(c::VN, r::N) where {N<:AbstractFloat,VN<:AbstractVector{N}}
    return Ball2(c, r)
end

include("init.jl")

end  # module
