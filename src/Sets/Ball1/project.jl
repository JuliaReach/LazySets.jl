function project(B::Ball1, block::AbstractVector{Int}; kwargs...)
    return Ball1(B.center[block], B.radius)
end
