@validate function project(B::BallInf, block::AbstractVector{Int}; kwargs...)
    return BallInf(B.center[block], B.radius)
end
