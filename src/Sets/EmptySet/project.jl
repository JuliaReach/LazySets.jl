function project(::EmptySet{N}, block::AbstractVector{Int}; kwargs...) where {N}
    return EmptySet{N}(length(block))
end
