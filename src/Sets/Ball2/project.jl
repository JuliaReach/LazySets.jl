@validate function project(B::Ball2, block::AbstractVector{Int}; kwargs...)
    return Ball2(B.center[block], B.radius)
end
