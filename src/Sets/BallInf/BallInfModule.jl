module BallInfModule

using Reexport

using ..LazySets: AbstractHyperrectangle, _high_AbstractBallp, _low_AbstractBallp,
                  _ρ_sev_hyperrectangle, _σ_sev_hyperrectangle
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Arrays: SingleEntryVector
using ReachabilityBase.Comparison: isapproxzero
using ReachabilityBase.Distribution: reseed!
using LinearAlgebra: dot

@reexport import ..API: area, center, isoperationtype, high, low, radius, rand,
                        reflect, volume, project, scale, ρ, σ, translate!
@reexport import ..LazySets: ball_norm, isflat, ngens, radius_ball,
                             radius_hyperrectangle, □
@reexport using ..API

export BallInf

include("BallInf.jl")

include("area.jl")
include("ball_norm.jl")
include("center.jl")
include("genmat.jl")
include("high.jl")
include("isflat.jl")
include("isoperationtype.jl")
include("low.jl")
include("ngens.jl")
include("radius.jl")
include("radius_ball.jl")
include("radius_hyperrectangle.jl")
include("rand.jl")
include("reflect.jl")
include("volume.jl")
include("project.jl")
include("scale.jl")
include("support_function.jl")
include("support_vector.jl")
include("translate.jl")

function □(c::VN, r::N) where {N,VN<:AbstractVector{N}}
    return BallInf(c, r)
end

end  # module