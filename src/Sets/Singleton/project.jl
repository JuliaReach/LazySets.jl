@validate function project(S::Singleton, block::AbstractVector{Int}; kwargs...)
    return Singleton(center(S)[block])
end
