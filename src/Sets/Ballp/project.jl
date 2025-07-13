@validate function project(B::Ballp, block::AbstractVector{Int}; kwargs...)
    return Ballp(B.p, B.center[block], B.radius)
end
