@validate function project(cup::UnionSetArray, block::AbstractVector{Int}; kwargs...)
    return UnionSetArray([project(X, block; kwargs...) for X in cup])
end
