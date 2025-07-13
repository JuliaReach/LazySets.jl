@validate function project(U::Universe{N}, block::AbstractVector{Int}; kwargs...) where {N}
    return Universe{N}(length(block))
end
