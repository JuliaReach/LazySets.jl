function project(∅::EmptySet{N}, block::AbstractVector{Int}; kwargs...) where {N}
    @assert length(block) <= dim(∅) "incompatible dimensions $(length(block)) and $(dim(∅))"
    @assert all(1 <= v <= dim(∅) for v in block) "invalid dimension in index vector"

    return EmptySet{N}(length(block))
end
