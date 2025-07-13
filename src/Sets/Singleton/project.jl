@validate function project(S::Singleton, block::AbstractVector{Int}; kwargs...)
    return Singleton(element(S)[block])
end
